/*****************************************************************************
 TVBii.c
 The Virtual Brain
 (c) 2017 - Michael Schirner
 BrainModes Group
 Charité University Medicine Berlin
 m.schirner@fu-berlin.de
 Licensed under the GNU General Public License 3.0 license.
 ******************************************************************************/

#include <stdio.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

struct Xi_p{
    float **Xi_elems;
};

struct SC_capS{
    float *cap;
};

FILE *FCout, *WFout, *Sout;

/*
 Open files for writing simulated activity
 */
void openFCoutfile(char *paramset, int myid, char *sub_id){
    char outfilename[1000];memset(outfilename, 0, 1000*sizeof(char));
    char buffer[10];memset(buffer, 0, 10*sizeof(char));
    char underscore[2];
    underscore[0]='_';
    underscore[1]='\0';
    strcpy (outfilename,"output/");
    strcat (outfilename,"/BOLD_");
    strcat (outfilename,sub_id);strcat (outfilename,underscore);
    strcat (outfilename,paramset);strcat (outfilename,underscore);
    snprintf(buffer, 10, "%d", myid);strcat (outfilename,buffer);
    strcat (outfilename,".txt");
    FCout = fopen(outfilename, "w");
    /*
    memset(outfilename, 0, 1000*sizeof(char));
    strcpy (outfilename,"output/");
    strcat (outfilename,"/INP_");
    strcat (outfilename,sub_id);strcat (outfilename,underscore);
    strcat (outfilename,paramset);strcat (outfilename,underscore);
    snprintf(buffer, 10, "%d", myid);strcat (outfilename,buffer);
    strcat (outfilename,".txt");
    WFout = fopen(outfilename, "w");
    memset(outfilename, 0, 1000*sizeof(char));
    strcpy (outfilename,"output/");
    strcat (outfilename,"/S_");
    strcat (outfilename,sub_id);strcat (outfilename,underscore);
    strcat (outfilename,paramset);strcat (outfilename,underscore);
    snprintf(buffer, 10, "%d", myid);strcat (outfilename,buffer);
    strcat (outfilename,".txt");
    Sout = fopen(outfilename, "w");
     */
}


/*
 Import SC (long-range coupling) and convert into efficient data structure
 */
int importGlobalConnectivity(char *SC_cap_filename, char *SC_dist_filename, char *SC_inputreg_filename, int regions, float **region_activity, struct Xi_p **reg_globinp_p, float global_trans_v, int **n_conn_table, float **n_conn_table_G_NMDA, float G_J_NMDA, struct SC_capS **SC_cap, float **SC_rowsums)
{
    
    int i,j,k, tmp3, maxdelay=0, tmpint;
    float *region_activity_p;
    double tmp, tmp2;
    struct Xi_p *reg_globinp_pp;
    struct SC_capS *SC_capp;
    int num_incoming_conn_cap=0, num_incoming_conn_dist=0, num_incoming_conn_inputreg=0;
    
    // Open SC files
    FILE *file_cap, *file_dist, *file_inputreg;
    file_cap=fopen(SC_cap_filename, "r");
    file_dist=fopen(SC_dist_filename, "r");
    file_inputreg=fopen(SC_inputreg_filename, "r");
    if (file_cap==NULL || file_dist==NULL || file_inputreg==NULL)
    {
        printf( "\n ERROR: Could not open SC files. Terminating... \n\n");
        exit(0);
    }
    
    // Read number of regions in header and check whether it's consistent with other specifications
    int readSC_cap, readSC_dist, readSC_inp;
    if(fscanf(file_cap,"%d",&readSC_cap) == EOF || fscanf(file_dist,"%d",&readSC_dist) == EOF || fscanf(file_inputreg,"%d",&readSC_inp) == EOF){
        printf( "\n ERROR: Unexpected end-of-file in file. File contains less input than expected. Terminating... \n\n");
        exit(0);
    }
    if (readSC_cap != regions || readSC_dist != regions || readSC_inp != regions) {
        printf( "\n ERROR: Inconsistent number of large-scale regions in SC files. Terminating... \n\n");
        fclose(file_dist);fclose(file_inputreg);
        exit(0);
    }
    
    // Allocate a counter (to store number of region inputs for each region) and the SCcap array (to store connection weights)
    *SC_rowsums = (float *)_mm_malloc(regions*sizeof(float),16);
    *n_conn_table = (int *)_mm_malloc(regions*sizeof(int),16);
    *n_conn_table_G_NMDA = (float *)_mm_malloc(regions*sizeof(float),16);
    *SC_cap = (struct SC_capS *)_mm_malloc(regions*sizeof(struct SC_capS),16);
    SC_capp = *SC_cap;
    if(*n_conn_table==NULL || *n_conn_table_G_NMDA==NULL || SC_capp==NULL || SC_rowsums==NULL){
        printf("Running out of memory. Terminating.\n");fclose(file_dist);fclose(file_cap);fclose(file_inputreg);exit(2);
    }
    
    // Read the maximal fiber length in header of SCdist-file and compute maxdelay
    if(fscanf(file_dist,"%lf",&tmp) == EOF){
        printf( "ERROR: Unexpected end-of-file in file. File contains less input than expected. Terminating... \n");
        exit(0);
    }
    maxdelay = (int)((tmp/global_trans_v)+0.5); // +0.5 for rounding by casting
    if (maxdelay < 1) maxdelay = 1; // Case: no time delays
    
    // Allocate ringbuffer that contains region activity for each past time-step until maxdelay and another ringbuffer that contains pointers to the first ringbuffer
    *region_activity = (float *)_mm_malloc(maxdelay*regions*sizeof(float),16);
    region_activity_p = *region_activity;
    *reg_globinp_p = (struct Xi_p *)_mm_malloc(maxdelay*regions*sizeof(struct Xi_p),16);
    reg_globinp_pp = *reg_globinp_p;
    if(region_activity_p==NULL || reg_globinp_p==NULL){
        printf("Running out of memory. Terminating.\n");fclose(file_dist);exit(2);
    }
    for (j=0; j<maxdelay*regions; j++) {
        region_activity_p[j]=0.001;
    }
    
    // Read SC files and set pointers for each input region and correspoding delay for each ringbuffer time-step
    int ring_buff_position;
    for (i=0; i<regions; i++) {
        // Read region index of current region (first number of each row) and check whether its consistent for all files
        if(fscanf(file_cap,"%d",&num_incoming_conn_cap) == EOF || fscanf(file_dist,"%d",&num_incoming_conn_dist) == EOF || fscanf(file_inputreg,"%d",&num_incoming_conn_inputreg) == EOF){
            printf( "ERROR: Unexpected end-of-file in SC files. File(s) contain(s) less input than expected. Terminating... \n");
            exit(0);
        }
        if (num_incoming_conn_cap != i || num_incoming_conn_dist != i || num_incoming_conn_inputreg != i) {
            printf( "ERROR: Inconsistencies in global input files, seems like row number is incorrect in some files. i=%d cap=%d dist=%d inp=%d Terminating.\n\n", i, num_incoming_conn_cap, num_incoming_conn_dist, num_incoming_conn_inputreg);
            exit(0);
        }
        
        // Read number of incoming connections for this region (second number of each row) and check whether its consistent across input files
        if(fscanf(file_cap,"%d",&num_incoming_conn_cap) == EOF || fscanf(file_dist,"%d",&num_incoming_conn_dist) == EOF || fscanf(file_inputreg,"%d",&num_incoming_conn_inputreg) == EOF){
            printf( "ERROR: Unexpected end-of-file in file %s or %s. File contains less input than expected. Terminating... \n\n", SC_dist_filename, SC_inputreg_filename);
            exit(0);
        }
        if (num_incoming_conn_cap != num_incoming_conn_inputreg || num_incoming_conn_dist != num_incoming_conn_inputreg) {
            printf( "ERROR: Inconsistencies in SC files: Different numbers of input connections. Terminating.\n\n");
            exit(0);
        }
        
        (*n_conn_table)[i]      = num_incoming_conn_inputreg;
        if (num_incoming_conn_inputreg > 0) {
            (*n_conn_table_G_NMDA)[i]     = G_J_NMDA  / num_incoming_conn_inputreg;
            //(*n_conn_table_G_NMDA)[i]     = G_J_NMDA;
        } else{
            (*n_conn_table_G_NMDA)[i]     = 0;
        }
        
        
        if ((*n_conn_table)[i] > 0) {
            // SC strength
            SC_capp[i].cap = (float *)_mm_malloc(((*n_conn_table)[i])*sizeof(float),16);
            if(SC_capp[i].cap==NULL){
                printf("Running out of memory. Terminating.\n");exit(2);
            }
            
            // Allocate memory for input-region-pointer arrays for each time-step in ringbuffer
            for (j=0; j<maxdelay; j++){
                reg_globinp_pp[i+j*regions].Xi_elems=(float **)_mm_malloc(((*n_conn_table)[i])*sizeof(float *),16);
                if(reg_globinp_pp[i+j*regions].Xi_elems==NULL){
                    printf("Running out of memory. Terminating.\n");exit(2);
                }
            }
            
            float sum_caps=0.0;
            // Read incoming connections and set pointers
            for (j=0; j<(*n_conn_table)[i]; j++) {
                
                if(fscanf(file_cap,"%lf",&tmp) != EOF && fscanf(file_dist,"%lf",&tmp2) != EOF && fscanf(file_inputreg,"%d",&tmp3) != EOF){
                    
                    tmpint = (int)((tmp2/global_trans_v)+0.5); // +0.5 for rounding by casting
                    if (tmpint < 0){
                        printf( "\nERROR: Negative connection length/delay %d -> %d. Terminating... \n\n",i,tmp3);exit(0);
                    }
                    
                    //SC_capp[i].cap[j] = (float)tmp * (*n_conn_table_G_NMDA)[i];
                    SC_capp[i].cap[j] = (float)tmp * G_J_NMDA;
                    sum_caps += (float)tmp;
                    
                    if (tmp3 >= 0) {
                        ring_buff_position=maxdelay*regions - tmpint*regions + tmp3;
                        for (k=0; k<maxdelay; k++) {
                            reg_globinp_pp[i+k*regions].Xi_elems[j]=&region_activity_p[ring_buff_position];
                            ring_buff_position += regions;
                            if (ring_buff_position > (maxdelay*regions-1)) ring_buff_position -= maxdelay*regions;
                        }
                    } else {
                        printf( "\nERROR: Region index is negative: %d -> %d. Terminating... \n\n",i,tmp3);exit(0);
                    }
                    
                    
                } else{
                    printf( "\n ERROR: Unexpected end-of-file in file %s or %s. File contains less input than expected. Terminating... \n\n", SC_inputreg_filename, SC_dist_filename);
                    exit(0);
                }
                
            }
            if (sum_caps <= 0) {
                printf( "\nERROR: Sum of connection strenghts is negative or zero. sum-caps node %d = %f. Terminating... \n\n",i,sum_caps);exit(0);
            }
            (*SC_rowsums)[i] = sum_caps;
        }
    }
    
    fclose(file_dist);fclose(file_inputreg);
    return maxdelay;
}




int main(int argc, char *argv[])
{
    
   
    /*************
     Set up model
     *************/
    
    int myid=0, i, j, k;
    
    //Get starting time
    time_t start = time(NULL);
    
    // Check input arguments
    if (argc != 4 || atoi(argv[3]) > 8  || atoi(argv[3]) < 0) {
        printf( "\n ERROR: Too few/many or wrong arguments. Terminating... \n\n");
        exit(0);
    }
    
    // Open output files and switch between different types of input activity
    char outputfolder[100];memset(outputfolder, 0, 100*sizeof(char));
    char subject_file[100];memset(subject_file, 0, 100*sizeof(char));
    strcpy(subject_file,"input/");strcat(subject_file,argv[2]);
    switch (atoi(argv[3])) {
        case 0:
            strcat(subject_file,"_fullextinp1000.txt");
            strcpy(outputfolder,"full");
            break;
        case 1:
            strcat(subject_file,"_fullextinprandperm1000.txt");
            strcpy(outputfolder,"fullrand");
            break;
        case 2:
            strcat(subject_file,"_fullextinp_diff1000.txt");
            strcpy(outputfolder,"diff");
            break;
        case 3:
            strcat(subject_file,"_fullextinprandperm_diff1000.txt");
            strcpy(outputfolder,"diffrand");
            break;
        case 4:
            strcat(subject_file,"_fullextinp1000.txt");
            strcpy(outputfolder,"testfull");
            break;
        case 5:
            strcat(subject_file,"_fullextinp_diff1000.txt");
            strcpy(outputfolder,"testdiff");
            break;
        case 6:
            memset(subject_file, 0, 100*sizeof(char));
            strcpy(subject_file,"input/");
            strcat(subject_file,"artifical_powermod_input.txt");
            strcpy(outputfolder,"testdiff");
            break;
        case 7:
            memset(subject_file, 0, 100*sizeof(char));
            strcpy(subject_file,"input/");
            strcat(subject_file,"artifical_powermod_input_9Hz.txt");
            strcpy(outputfolder,"testdiff");
            break;
        case 8:
            memset(subject_file, 0, 100*sizeof(char));
            strcpy(subject_file,"input/");
            strcat(subject_file,"artifical_powermod_input_11Hz.txt");
            strcpy(outputfolder,"testdiff");
            break;
    }
    openFCoutfile(argv[1], myid, argv[2]);
    
    
    
    /*
     Global model and integration parameters
     */
    float TR                      = 1.94;             // (s)  TR of fMRI data
    int   TR_ts                   = TR * 1000;        // n-th time step for writing out a fMRI data point
    const float dt                = 0.1;              //      Integration step length dt = 0.1 ms
    const float model_dt          = 0.001;            // (Hz) Time-step of model (sampling-rate=1000 Hz)
    const int vectorization_grade = 4;                // depends on CPU Architecture and available intrinsics.
    int   time_steps              = 667*TR*1000;      // simulation length: 667 TRs * 1.94 s/TR * 1000 time-steps/s
    int   time_steps_FIC_burnin   = 10*1000;          // 10 s as a burn-in for FIC tuning to let model activity stabilized
    int   time_steps_FIC_fit      = time_steps - time_steps_FIC_burnin;
    const int   nodes             = 68;               // Number of nodes; must be a multiple of vectorization grade
    const int   nodes_vec         = nodes/vectorization_grade;
    const int   regions           = 68;               // Number of large-scale regions (identical to nodes for non-surface simulations)
    float global_trans_v          = 1.0;              // Global transmission velocity (m/s)
    float G                       = 0.5;              // Global coupling strength
    int   FIC_iters               = 12;               // Number of FIC tuning iterations
    const float avg_FR            = 3.0631;           // (Hz) avg firing rate of an exc. population
    const float tuning_factor     = 0.005;            // parameter for FIC tuning
    float adaptive_tuning_factor  = tuning_factor;    // parameter for FIC tuning
    const float meanFRfact        = 1.0 / (time_steps_FIC_fit*10); // Factor for averaging firing rates over FIC simulation period
    float mean_mean_FR            = 0.0;              // variable for FIC tuning
    float variance_FR             = 0.0;              // variable for FIC tuning
    float tmpglobinput, tmpBOLD;
    int   ring_buf_pos=0;                             // Position in ring buffer of past state variables for long-range input
    float tmp_exp_E[4]          __attribute__((aligned(16)));
    float tmp_exp_I[4]          __attribute__((aligned(16)));
    float tmp_E[4]          __attribute__((aligned(16)));
    float tmp_I[4]          __attribute__((aligned(16)));
    
    
    /*
     Local model: DMF-Parameters from Deco et al. JNeuro 2014
     */
    float w_plus  = 1.4;                // local excitatory recurrence
    //float I_ext   = 0;                // Additional external stimulation -- not used in hybrid models
    //float J_NMDA  = 0.15;             // (nA) excitatory synaptic coupling
    float J_NMDA  = 0.0;                // J_NMDA is replaced by injected input in hybrid models
    //float J_i     = 1.0;              // Single value for no-FIC models
    const float a_E     = 310;          // (n/C)
    const float b_E     = 125;          // (Hz)
    const float d_E     = 0.16;         // (s)
    const float a_I     = 615;          // (n/C)
    const float b_I     = 177;          // (Hz)
    const float d_I     = 0.087;        // (s)
    const float gamma   = 0.641/1000.0; // factor 1000 for expressing everything in ms
    const float tau_E   = 100;          // (ms) Time constant of NMDA (excitatory)
    const float tau_I   = 10;           // (ms) Time constant of GABA (inhibitory)
    const float I_0     = 0.382;        // (nA) overall effective external input
    const float w_E     = 1.0;          // scaling of external input for excitatory pool
    const float w_I     = 0.7;          // scaling of external input for inhibitory pool
    const float gamma_I = 1.0/1000.0;   // for expressing inhib. pop. in ms
    float       eIf_e     = 1.0;        // weighting factor for external input (exc. pop.)
    float       eIf_i     = 1.0;        // weighting factor for external input (inh. pop.)
    // Derived parameters for more efficient computation
    const float min_d_E       = -1.0 * d_E;
    const float min_d_I       = -1.0 * d_I;
    const float imintau_E     = -1.0 / tau_E;
    const float imintau_I     = -1.0 / tau_I;
    const float w_E__I_0      = w_E * I_0;
    const float w_I__I_0      = w_I * I_0;
    const float one           = 1.0;
    const float w_plus__J_NMDA= w_plus * J_NMDA;
    
    
    /*
     Allocate and Initialize arrays
     */
    // Brain network model arrays
    float *ext_input = (float *)_mm_malloc(time_steps*nodes * sizeof(float),16);
    double tmp;
    float *S_i_E                    = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *S_i_I                    = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *r_i_E                    = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *r_i_I                    = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *EP_selfExc               = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *EP_globInp               = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *EP_inhInp                = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *EP_extInp                = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *EP_extInpinh             = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *global_input             = (float *)_mm_malloc(nodes * sizeof(float),16);
    float *J_i                      = (float *)_mm_malloc(nodes * sizeof(float),16);  // (nA) inhibitory synaptic coupling
    float *meanFR                   = (float *)_mm_malloc(nodes * sizeof(float),16);  // summation array for mean firing rate
    //Balloon-Windkessel model arrays (see Deco et al. 2014 JNeuro)
    int   output_ts = time_steps / (TR / model_dt);                             // Length of BOLD time-series written to HDD
    int   num_output_ts      = 68;                                              // Number of BOLD time-series that are writte to HDD
    float *bw_x_ex    = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 1 of BW-model (exc. pop.)
    float *bw_f_ex    = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 2 of BW-model (exc. pop.)
    float *bw_nu_ex   = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 3 of BW-model (exc. pop.)
    float *bw_q_ex    = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 4 of BW-model (exc. pop.)
    float *bw_x_in    = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 1 of BW-model (inh. pop.)
    float *bw_f_in    = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 2 of BW-model (inh. pop.)
    float *bw_nu_in   = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 3 of BW-model (inh. pop.)
    float *bw_q_in    = (float *)_mm_malloc(num_output_ts * sizeof(float),16);  // State-variable 4 of BW-model (inh. pop.)
    if (S_i_E == NULL || S_i_I == NULL || r_i_E == NULL || r_i_I == NULL ||  global_input == NULL || J_i == NULL || meanFR == NULL || EP_selfExc == NULL || EP_globInp == NULL || EP_inhInp == NULL || EP_extInp == NULL || EP_extInpinh == NULL || ext_input == NULL|| bw_x_ex == NULL|| bw_f_ex == NULL|| bw_nu_ex == NULL|| bw_q_ex == NULL|| bw_q_in == NULL) {
        printf( "ERROR: Running out of memory. Aborting... \n");
        _mm_free(S_i_E);_mm_free(S_i_I);_mm_free(global_input);_mm_free(J_i);_mm_free(meanFR);
        return 1;
    }
    
    // Initialize model state variables / parameters
    for (j = 0; j < nodes; j++) {
        S_i_E[j]            = 0.001;
        S_i_I[j]            = 0.001;
        global_input[j]     = 0.001;
        meanFR[j]           = 0.0;
        J_i[j]              = 1.0;
    }
    
    // Initialize Balloon-Windkessel model parameters and arrays
    for (j = 0; j < num_output_ts; j++) {
        bw_x_ex[j] = 0.0;
        bw_f_ex[j] = 1.0;
        bw_nu_ex[j] = 1.0;
        bw_q_ex[j] = 1.0;
        bw_x_in[j] = 0.0;
        bw_f_in[j] = 1.0;
        bw_nu_in[j] = 1.0;
        bw_q_in[j] = 1.0;
    }
    float BOLD_ex[num_output_ts][output_ts+2];
    float rho = 0.34, alpha = 0.32, tau = 0.98, y = 1.0/0.41, kappa = 1.0/0.65;
    float V_0 = 0.02, k1 = 7 * rho, k2 = 2.0, k3 = 2 * rho - 0.2, ialpha = 1.0/alpha, itau = 1.0/tau, oneminrho = (1.0 - rho);
    float f_tmp;
    int   BOLD_offset = 11;
    int   BOLD_len_i  = -1;
    
    
    /*
     Load external input -- to be injected into hybrid model
     */
    FILE *ext_input_file = fopen(subject_file,"r");
    if (ext_input_file==NULL)    {
        printf( "\n ERROR: Could not open external input file. Terminating... \n\n");
        exit(0);
    }
    for (i = 0 ; i < time_steps*nodes; i++){
        if(fscanf(ext_input_file,"%lf",&tmp) != EOF){
            ext_input[i] = (float)tmp;
        } else{
            printf( "\n ERROR: Unexpected end-of-file in external input file. File contains less input than expected. Terminating... \n\n");
            exit(0);
        }
    }
    fclose(ext_input_file);
    
    
    
    /*
     Input-file mode is on: override some parameter-values as specified in additional param-file
     */
    float tmpratio = 0.0;
    if (argc > 1) {
        FILE *file;
        file=fopen(argv[1], "r");
        if (file==NULL)
        {
            printf( "\n ERROR: Could not open file %s. Terminating... \n\n", argv[1]);
            exit(0);
        }
        int skip;
        for (skip=-1; skip<myid; skip++) {
            if(fscanf(file,"%f",&G) != EOF && fscanf(file,"%f",&eIf_i) != EOF && fscanf(file,"%f",&tmpratio) != EOF){
            } else{
                printf( "\n ERROR: Unexpected end-of-file in file %s. File contains less input than expected. Terminating... \n\n", argv[1]);
                exit(0);
            }
        }
        fclose(file);
    }
    eIf_e = eIf_i / tmpratio;
    //const float G_J_NMDA      = G * J_NMDA;
    const float G_J_NMDA      = G ; // J_NMDA is replaced by injected input in hybrid models
    
    
    /*
     Import and set-up long-range connectivity
     */
    int         *n_conn_table;          // Number of incoming connections for each node
    float       *region_activity;       // Ring-buffer of past values of S_E for computing long-range input
    float       *n_conn_table_G_NMDA;   // legacy
    float       *SC_rowsums;            // legacy
    struct Xi_p *reg_globinp_p;         // Array of addresses poiting into region_activity buffer
    struct SC_capS *SC_cap;             // structural connectivity weights
    
    // Construct filenames
    char cap_file[100];memset(cap_file, 0, 100*sizeof(char));
    strcpy(cap_file,"input/");strcat(cap_file,argv[2]);strcat(cap_file,"_SCcap_large.txt");
    char dist_file[100];memset(dist_file, 0, 100*sizeof(char));
    strcpy(dist_file,"input/");strcat(dist_file,argv[2]);strcat(dist_file,"_SCdist_large.txt");
    char reg_file[100];memset(reg_file, 0, 100*sizeof(char));
    strcpy(reg_file,"input/");strcat(reg_file,argv[2]);strcat(reg_file,"_SCinputreg_large.txt");
    
    // Actual function to read and generate SC
    int         maxdelay = importGlobalConnectivity(cap_file, dist_file, reg_file, regions, &region_activity, &reg_globinp_p, global_trans_v, &n_conn_table, &n_conn_table_G_NMDA, G_J_NMDA, &SC_cap, &SC_rowsums);
    int         reg_act_size = regions * maxdelay; // size of region_activity ringbuffer
    
    
    /*
     Initialize or cast to vector-intrinsics types for required variables & parameters
     */
    const __m128    _dt                 = _mm_load1_ps(&dt);
    const __m128    _w_plus_J_NMDA      = _mm_load1_ps(&w_plus__J_NMDA);
    const __m128    _a_E                = _mm_load1_ps(&a_E);
    const __m128    _b_E                = _mm_load1_ps(&b_E);
    const __m128    _min_d_E            = _mm_load1_ps(&min_d_E);
    const __m128    _a_I                = _mm_load1_ps(&a_I);
    const __m128    _b_I                = _mm_load1_ps(&b_I);
    const __m128    _min_d_I            = _mm_load1_ps(&min_d_I);
    const __m128    _gamma              = _mm_load1_ps(&gamma);
    const __m128    _gamma_I            = _mm_load1_ps(&gamma_I);
    const __m128    _imintau_E          = _mm_load1_ps(&imintau_E);
    const __m128    _imintau_I          = _mm_load1_ps(&imintau_I);
    const __m128    _w_E__I_0           = _mm_load1_ps(&w_E__I_0);
    const __m128    _w_I__I_0           = _mm_load1_ps(&w_I__I_0);
    //const __m128    _I_0                = _mm_load1_ps(&I_0);
    const __m128    _one                = _mm_load1_ps(&one);
    const __m128    _J_NMDA             = _mm_load1_ps(&J_NMDA);
    const __m128    _eIf_e              = _mm_load1_ps(&eIf_e);
    const __m128    _eIf_i              = _mm_load1_ps(&eIf_i);
    float zero = 0.0, pone = 0.1;
    const __m128    _zero               = _mm_load1_ps(&zero);
    const __m128    _pone               = _mm_load1_ps(&pone);
    __m128          *_S_i_E             = (__m128*)S_i_E;
    __m128          *_S_i_I             = (__m128*)S_i_I;
    __m128          *_r_i_E             = (__m128*)r_i_E;
    __m128          *_r_i_I             = (__m128*)r_i_I;
    __m128          *_EP_selfExc        = (__m128*)EP_selfExc;
    __m128          *_EP_globInp        = (__m128*)EP_globInp;
    __m128          *_EP_inhInp         = (__m128*)EP_inhInp;
    __m128          *_EP_extInp         = (__m128*)EP_extInp;
    __m128          *_EP_extInpinh      = (__m128*)EP_extInpinh;
    __m128          *_tmp_E             = (__m128*)tmp_E;
    __m128          *_tmp_I             = (__m128*)tmp_I;
    __m128          *_tmp_exp_E         = (__m128*)tmp_exp_E;
    __m128          *_tmp_exp_I         = (__m128*)tmp_exp_I;
    __m128          *_global_input      = (__m128*)global_input;
    __m128          *_ext_input         = (__m128*)ext_input;
    __m128          *_J_i               = (__m128*)J_i;
    __m128          *_meanFR            = (__m128*)meanFR;
    __m128          _tmp_I_E, _tmp_I_I;
    __m128          _tmp_H_E, _tmp_H_I;
    
    
    
    
    
    /************
     FIC tuning
     ************/
    
    
    // Parameters for FIC tuning
    float previous_meanFR[nodes], Ji_save_bestFRdistance[nodes], meanFR_save_bestFRdistance[nodes];
    float avg_previous_Ji_change=0, avg_previous_FR_change, current_best_FRdistance=9999999.0, current_best_FRmean = 0.0,       previous_mean_meanFR=0,previous_adaptive_tuning_factor = tuning_factor, current_best_FR_SD=0.0, maximum_tuning_factor=-1;
    int is_improved = 0, is_fine_tuning=-1;
    int ts, int_i, i_node_vec, FIC_iter, ext_inp_counter=0;
    

    
    
    /*
     FIC tuning loop
     */
    for (FIC_iter = 0; FIC_iter < FIC_iters; FIC_iter++) {
        
        // Reset arrays
        ext_inp_counter=0;
        for (j = 0; j < nodes; j++) {
            S_i_E[j]            = 0.001;
            S_i_I[j]            = 0.001;
            global_input[j]     = 0.001;
            meanFR[j]           = 0.0;
        }
        /*
        for (j=0; j<reg_act_size; j++) {
            region_activity[j]=0.001;
        }
        */
        /*
        for (j = 0; j < nodes; j++){
            fprintf(FCout, "%.4f ",J_i[j]);
        }
        fprintf(FCout, "\n");
        */
        
        /*
         Phase I: Generate burn-in activity
         */
        for (ts = 0; ts < time_steps_FIC_burnin; ts++) {
            for (int_i = 0; int_i < 10; int_i++) {
                
                /*
                 Compute long-range coupling
                 */
                for(j=0; j<regions; j++){
                    tmpglobinput = 0;
                    for (k=0; k<n_conn_table[j]; k++) {
                        tmpglobinput += *reg_globinp_p[j+ring_buf_pos].Xi_elems[k] * SC_cap[j].cap[k];
                    }
                    global_input[j] = tmpglobinput;
                }
                
                /*
                 Integrate model
                 */
                for (i_node_vec = 0; i_node_vec < nodes_vec; i_node_vec++) {
                    // Excitatory population firing rate
                    _tmp_I_E    = _mm_sub_ps(_mm_mul_ps(_a_E,_mm_add_ps(_mm_add_ps(_mm_add_ps(_w_E__I_0,_mm_mul_ps(_eIf_e, _ext_input[ext_inp_counter])),_mm_mul_ps(_w_plus_J_NMDA, _S_i_E[i_node_vec])),_mm_sub_ps(_global_input[i_node_vec],_mm_mul_ps(_J_i[i_node_vec], _S_i_I[i_node_vec])))),_b_E);

                    *_tmp_exp_E   = _mm_mul_ps(_min_d_E, _tmp_I_E);

                    
                    tmp_exp_E[0]  = tmp_exp_E[0] != 0 ? expf(tmp_exp_E[0]) : 0.9; // To avoid division-by-zero
                    tmp_exp_E[1]  = tmp_exp_E[1] != 0 ? expf(tmp_exp_E[1]) : 0.9; // which can occur in
                    tmp_exp_E[2]  = tmp_exp_E[2] != 0 ? expf(tmp_exp_E[2]) : 0.9; // Eqs. (7+8) of
                    tmp_exp_E[3]  = tmp_exp_E[3] != 0 ? expf(tmp_exp_E[3]) : 0.9; // Deco et al. 2014 JNeuro
                    *_tmp_E       = _mm_div_ps(_tmp_I_E, _mm_sub_ps(_one, *_tmp_exp_E));

    
                    
                    _tmp_H_E  = *_tmp_E;
   
                    // Inhibitory population firing rate
                    _tmp_I_I = _mm_sub_ps(_mm_mul_ps(_a_I,_mm_sub_ps(_mm_add_ps(_mm_add_ps(_w_I__I_0,_mm_mul_ps(_eIf_i, _ext_input[ext_inp_counter])),_mm_mul_ps(_J_NMDA, _S_i_E[i_node_vec])), _S_i_I[i_node_vec])),_b_I);
                    ext_inp_counter++;
                    *_tmp_exp_I   = _mm_mul_ps(_min_d_I, _tmp_I_I);
                    tmp_exp_I[0]  = tmp_exp_I[0] != 0 ? expf(tmp_exp_I[0]) : 0.9;
                    tmp_exp_I[1]  = tmp_exp_I[1] != 0 ? expf(tmp_exp_I[1]) : 0.9;
                    tmp_exp_I[2]  = tmp_exp_I[2] != 0 ? expf(tmp_exp_I[2]) : 0.9;
                    tmp_exp_I[3]  = tmp_exp_I[3] != 0 ? expf(tmp_exp_I[3]) : 0.9;
                    *_tmp_I  = _mm_div_ps(_tmp_I_I, _mm_sub_ps(_one, *_tmp_exp_I));
                    
                    
                    _tmp_H_I  = *_tmp_I;
                    
                    _r_i_E[i_node_vec] =  _tmp_H_E;
                    _r_i_I[i_node_vec] =  _tmp_H_I;
                    
                    // Synaptic activity
                    _S_i_I[i_node_vec] = _mm_add_ps(_S_i_I[i_node_vec],_mm_mul_ps(_dt,_mm_add_ps(_mm_mul_ps(_imintau_I, _S_i_I[i_node_vec]),_mm_mul_ps(_tmp_H_I,_gamma_I))));
                    _S_i_E[i_node_vec] = _mm_add_ps(_S_i_E[i_node_vec],_mm_mul_ps(_dt, _mm_add_ps(_mm_mul_ps(_imintau_E, _S_i_E[i_node_vec]),_mm_mul_ps(_mm_mul_ps(_mm_sub_ps(_one, _S_i_E[i_node_vec]),_gamma),_tmp_H_E))));
                }
                
                for(j=0; j<nodes; j++){
                    S_i_E[j] = S_i_E[j] < 1 ? S_i_E[j] : 1;
                    S_i_E[j] = S_i_E[j] > 0 ? S_i_E[j] : 0;
                    S_i_I[j] = S_i_I[j] < 1 ? S_i_I[j] : 1;
                    S_i_I[j] = S_i_I[j] > 0 ? S_i_I[j] : 0;
                }
                
                
                ext_inp_counter -= nodes_vec;
                memcpy(&region_activity[ring_buf_pos], S_i_E, regions*sizeof( float ));
                // Shift ring-buff-pos
                ring_buf_pos = ring_buf_pos<(reg_act_size-regions) ? (ring_buf_pos+regions) : 0;
                
                /*
                if (FIC_iter == 4) {
                    for (j = 0; j < nodes; j++){
                        fprintf(WFout, "%.3f ",r_i_E[j]);
                    }
                    for (j = 0; j < nodes; j++){
                        fprintf(WFout, "%.3f ",r_i_I[j]);
                    }
                    fprintf(WFout, "\n");
                }
                 */
            }
            ext_inp_counter += nodes_vec;
        }
        /*
        if (FIC_iter == 4) {
            exit(0);
        }
        */
        /*
         Phase II: Simulate full time series and compute mean firing rate for each node
         */
        for (; ts < time_steps; ts++) {
            for (int_i = 0; int_i < 10; int_i++) {
                
                /*
                 Compute long-range coupling
                 */
                for(j=0; j<regions; j++){
                    tmpglobinput = 0;
                    for (k=0; k<n_conn_table[j]; k++) {
                        tmpglobinput += *reg_globinp_p[j+ring_buf_pos].Xi_elems[k] * SC_cap[j].cap[k];
                    }
                    global_input[j] = tmpglobinput;
                }
                
                /*
                 Integrate model
                 */
                for (i_node_vec = 0; i_node_vec < nodes_vec; i_node_vec++) {
                    // Excitatory population firing rate
                    _tmp_I_E    = _mm_sub_ps(_mm_mul_ps(_a_E,_mm_add_ps(_mm_add_ps(_mm_add_ps(_w_E__I_0,_mm_mul_ps(_eIf_e, _ext_input[ext_inp_counter])),_mm_mul_ps(_w_plus_J_NMDA, _S_i_E[i_node_vec])),_mm_sub_ps(_global_input[i_node_vec],_mm_mul_ps(_J_i[i_node_vec], _S_i_I[i_node_vec])))),_b_E);
                    
                    *_tmp_exp_E   = _mm_mul_ps(_min_d_E, _tmp_I_E);
                    
                    
                    tmp_exp_E[0]  = tmp_exp_E[0] != 0 ? expf(tmp_exp_E[0]) : 0.9; // To avoid division-by-zero
                    tmp_exp_E[1]  = tmp_exp_E[1] != 0 ? expf(tmp_exp_E[1]) : 0.9; // which can occur in
                    tmp_exp_E[2]  = tmp_exp_E[2] != 0 ? expf(tmp_exp_E[2]) : 0.9; // Eqs. (7+8) of
                    tmp_exp_E[3]  = tmp_exp_E[3] != 0 ? expf(tmp_exp_E[3]) : 0.9; // Deco et al. 2014 JNeuro
                    *_tmp_E       = _mm_div_ps(_tmp_I_E, _mm_sub_ps(_one, *_tmp_exp_E));
                    
                    _tmp_H_E  = *_tmp_E;
                    _meanFR[i_node_vec] = _mm_add_ps(_meanFR[i_node_vec],_tmp_H_E);
                    
                    // Inhibitory population firing rate
                    _tmp_I_I = _mm_sub_ps(_mm_mul_ps(_a_I,_mm_sub_ps(_mm_add_ps(_mm_add_ps(_w_I__I_0,_mm_mul_ps(_eIf_i, _ext_input[ext_inp_counter])),_mm_mul_ps(_J_NMDA, _S_i_E[i_node_vec])), _S_i_I[i_node_vec])),_b_I);
                    ext_inp_counter++;
                    *_tmp_exp_I   = _mm_mul_ps(_min_d_I, _tmp_I_I);
                    tmp_exp_I[0]  = tmp_exp_I[0] != 0 ? expf(tmp_exp_I[0]) : 0.9;
                    tmp_exp_I[1]  = tmp_exp_I[1] != 0 ? expf(tmp_exp_I[1]) : 0.9;
                    tmp_exp_I[2]  = tmp_exp_I[2] != 0 ? expf(tmp_exp_I[2]) : 0.9;
                    tmp_exp_I[3]  = tmp_exp_I[3] != 0 ? expf(tmp_exp_I[3]) : 0.9;
                    *_tmp_I  = _mm_div_ps(_tmp_I_I, _mm_sub_ps(_one, *_tmp_exp_I));
                    
                    _tmp_H_I  = *_tmp_I;
                    
                    _r_i_E[i_node_vec] =  _tmp_H_E;
                    _r_i_I[i_node_vec] =  _tmp_H_I;
                    
                    // Synaptic activity
                    _S_i_I[i_node_vec] = _mm_add_ps(_S_i_I[i_node_vec],_mm_mul_ps(_dt,_mm_add_ps(_mm_mul_ps(_imintau_I, _S_i_I[i_node_vec]),_mm_mul_ps(_tmp_H_I,_gamma_I))));
                    _S_i_E[i_node_vec] = _mm_add_ps(_S_i_E[i_node_vec],_mm_mul_ps(_dt, _mm_add_ps(_mm_mul_ps(_imintau_E, _S_i_E[i_node_vec]),_mm_mul_ps(_mm_mul_ps(_mm_sub_ps(_one, _S_i_E[i_node_vec]),_gamma),_tmp_H_E))));
                }
                for(j=0; j<nodes; j++){
                    S_i_E[j] = S_i_E[j] < 1 ? S_i_E[j] : 1;
                    S_i_E[j] = S_i_E[j] > 0 ? S_i_E[j] : 0;
                    S_i_I[j] = S_i_I[j] < 1 ? S_i_I[j] : 1;
                    S_i_I[j] = S_i_I[j] > 0 ? S_i_I[j] : 0;
                }
                
                ext_inp_counter -= nodes_vec;
                memcpy(&region_activity[ring_buf_pos], S_i_E, regions*sizeof( float ));
                // Shift ring-buff-pos
                ring_buf_pos = ring_buf_pos<(reg_act_size-regions) ? (ring_buf_pos+regions) : 0;
                

            }
            ext_inp_counter += nodes_vec;
        }

        
        /*
         Phase III: Evaluate mean FR and adapt J_i
         */
        float tmpFR=0.0;
        mean_mean_FR = 0.0;
        avg_previous_FR_change = 0.0;
        for (j = 0; j < nodes; j++){
            tmpFR = meanFR[j]*meanFRfact;
            meanFR[j] = tmpFR;
            mean_mean_FR += meanFR[j];
            avg_previous_FR_change += fabsf(tmpFR - previous_meanFR[j]);
            previous_meanFR[j] = tmpFR;
        }
        avg_previous_FR_change /= nodes;
        mean_mean_FR /= nodes;
        
        printf("%d %f\n",FIC_iter, mean_mean_FR);
        
        variance_FR = 0.0;
        for (j = 0; j < nodes; j++){
            variance_FR += (meanFR[j] - mean_mean_FR)*(meanFR[j] - mean_mean_FR);
        }
        variance_FR /= (nodes - 1);
        
        if (fabsf(mean_mean_FR - avg_FR) < current_best_FRdistance && maximum_tuning_factor < 0){
            previous_adaptive_tuning_factor = adaptive_tuning_factor;
        }
        if (FIC_iter > 0 && is_improved == 0 && maximum_tuning_factor < 0) {
            if (avg_previous_Ji_change != 0 && (previous_mean_meanFR - mean_mean_FR) != 0) {
                adaptive_tuning_factor = fabsf(avg_previous_Ji_change / (previous_mean_meanFR - mean_mean_FR));
            } else {
                adaptive_tuning_factor = tuning_factor;
            }
        } else if (FIC_iter == 0 && mean_mean_FR <= 2.7){
            for (j = 0; j < nodes; j++){
                J_i[j] = 0.5;
            }
        }
        
        if ((fabsf(mean_mean_FR - avg_FR) >= current_best_FRdistance || mean_mean_FR <= 2.7) && is_fine_tuning == -1){
            if (mean_mean_FR <= 2.7 && current_best_FRmean > avg_FR) {
                is_improved = 0;
                adaptive_tuning_factor = previous_adaptive_tuning_factor / 2;
                previous_adaptive_tuning_factor /= 2;
                maximum_tuning_factor = previous_adaptive_tuning_factor;
            } else{
                is_improved++;
            }
            
            for (j = 0; j < nodes; j++){
                J_i[j] = Ji_save_bestFRdistance[j];
                meanFR[j] = meanFR_save_bestFRdistance[j];
                previous_meanFR[j] = meanFR_save_bestFRdistance[j];
            }
        } else{
            is_improved = 0;
            if ((is_fine_tuning == 0 && current_best_FR_SD > sqrtf(variance_FR) && fabsf(mean_mean_FR - avg_FR) < 0.3) || is_fine_tuning == -1) {
                current_best_FRmean = mean_mean_FR;
                current_best_FRdistance = fabsf(mean_mean_FR - avg_FR);
                current_best_FR_SD = sqrtf(variance_FR);
                for (j = 0; j < nodes; j++){
                    Ji_save_bestFRdistance[j] = J_i[j];
                    meanFR_save_bestFRdistance[j] = meanFR[j];
                    previous_meanFR[j] = meanFR[j];
                }
                previous_mean_meanFR = mean_mean_FR;
            } else{
                is_fine_tuning = -1;
                for (j = 0; j < nodes; j++){
                    J_i[j] = Ji_save_bestFRdistance[j];
                    meanFR[j] = meanFR_save_bestFRdistance[j];
                    previous_meanFR[j] = meanFR_save_bestFRdistance[j];
                }
                is_improved++;
            }
            
        }
        
        if (adaptive_tuning_factor > 3*previous_adaptive_tuning_factor) {
            adaptive_tuning_factor = previous_adaptive_tuning_factor*2;
        } else if (adaptive_tuning_factor < previous_adaptive_tuning_factor/3){
            adaptive_tuning_factor = previous_adaptive_tuning_factor/2;
        }
        
        avg_previous_Ji_change = 0.0;
        float tmpJichange=0.0;
        if (is_improved == 0 && current_best_FRdistance >= 0.3) {
            for (j = 0; j < nodes; j++){
                tmpJichange = (meanFR[j] - avg_FR)*adaptive_tuning_factor;
//                if (J_i[j] + tmpJichange > 3) {
//                    tmpJichange = 3 - J_i[j];
//                }
                J_i[j] += tmpJichange;
                avg_previous_Ji_change += (tmpJichange);
                //printf("%.4f %.4f %.4f %.4f %.4f %.4f \n", meanFR[j], J_i[j], SC_rowsums[j], G, J_NMDA, eIf_e);
            }
            avg_previous_Ji_change /= nodes;
            is_fine_tuning = -1;
        } else{
            if (current_best_FRdistance < 0.3){
                is_fine_tuning = 0;
                for (j = 0; j < nodes; j++){
                    if (fabsf(meanFR_save_bestFRdistance[j] - avg_FR) > current_best_FR_SD*0.6) J_i[j] += (meanFR_save_bestFRdistance[j] - avg_FR)*adaptive_tuning_factor;
                }
            } else{
                is_fine_tuning = -1;
                if (current_best_FRmean < avg_FR) {
                    for (j = 0; j < nodes; j++){
                        if (meanFR_save_bestFRdistance[j] < (current_best_FRmean - current_best_FR_SD*0.6)) J_i[j] += (meanFR_save_bestFRdistance[j] - current_best_FRmean)*adaptive_tuning_factor;
                    }
                } else{
                    for (j = 0; j < nodes; j++){
                        if (meanFR_save_bestFRdistance[j] > (current_best_FRmean + current_best_FR_SD*0.6)) J_i[j] += (meanFR_save_bestFRdistance[j] - current_best_FRmean)*adaptive_tuning_factor;
                    }
                }
            }
        }
    } // end: FIC optimization
    
    // Set J_i to the best values that were found during optimization
    for (j = 0; j < nodes; j++){
        J_i[j] = Ji_save_bestFRdistance[j];
    }
    
    
    
    
    /************************************************
     The actual simulation with fitted J_i parameters
     *************************************************/
    int ts_bold=0;
    ext_inp_counter=0;
    for (ts = 0; ts < time_steps; ts++) {
        
        // Reset arrays for writing out state variables
        for (i_node_vec = 0; i_node_vec < nodes_vec; i_node_vec++) {
            _EP_selfExc[i_node_vec]     =_zero;
            _EP_globInp[i_node_vec]     =_zero;
            _EP_inhInp[i_node_vec]      =_zero;
            _EP_extInp[i_node_vec]      =_zero;
            _EP_extInpinh[i_node_vec]   =_zero;
        }
        
        for (int_i = 0; int_i < 10; int_i++) {
            /*
             Compute long-range coupling
             */
            for(j=0; j<regions; j++){
                tmpglobinput = 0;
                for (k=0; k<n_conn_table[j]; k++) {
                    tmpglobinput += *reg_globinp_p[j+ring_buf_pos].Xi_elems[k] * SC_cap[j].cap[k];
                }
                global_input[j] = tmpglobinput;
            }
            
            /*
             Integrate model
             */
            for (i_node_vec = 0; i_node_vec < nodes_vec; i_node_vec++) {
                
                // Store some state variables
                _EP_selfExc[i_node_vec]=_mm_add_ps(_EP_selfExc[i_node_vec],  _mm_mul_ps(_w_plus_J_NMDA, _S_i_E[i_node_vec]));
                _EP_globInp[i_node_vec]=_mm_add_ps(_EP_globInp[i_node_vec],  _global_input[i_node_vec]);
                _EP_inhInp[i_node_vec] =_mm_add_ps(_EP_inhInp[i_node_vec],   _mm_mul_ps(_J_i[i_node_vec], _S_i_I[i_node_vec]));
                _EP_extInp[i_node_vec] =_mm_add_ps(_EP_extInp[i_node_vec],   _mm_mul_ps(_eIf_e, _ext_input[ext_inp_counter]));
                _EP_extInpinh[i_node_vec] =_mm_add_ps(_EP_extInpinh[i_node_vec],   _mm_mul_ps(_eIf_i, _ext_input[ext_inp_counter]));
                
                // Excitatory population firing rate
                _tmp_I_E    = _mm_sub_ps(_mm_mul_ps(_a_E,_mm_add_ps(_mm_add_ps(_mm_add_ps(_w_E__I_0,_mm_mul_ps(_eIf_e, _ext_input[ext_inp_counter])),_mm_mul_ps(_w_plus_J_NMDA, _S_i_E[i_node_vec])),_mm_sub_ps(_global_input[i_node_vec],_mm_mul_ps(_J_i[i_node_vec], _S_i_I[i_node_vec])))),_b_E);
                
                *_tmp_exp_E   = _mm_mul_ps(_min_d_E, _tmp_I_E);
                
                
                tmp_exp_E[0]  = tmp_exp_E[0] != 0 ? expf(tmp_exp_E[0]) : 0.9; // To avoid division-by-zero
                tmp_exp_E[1]  = tmp_exp_E[1] != 0 ? expf(tmp_exp_E[1]) : 0.9; // which can occur in
                tmp_exp_E[2]  = tmp_exp_E[2] != 0 ? expf(tmp_exp_E[2]) : 0.9; // Eqs. (7+8) of
                tmp_exp_E[3]  = tmp_exp_E[3] != 0 ? expf(tmp_exp_E[3]) : 0.9; // Deco et al. 2014 JNeuro
                *_tmp_E       = _mm_div_ps(_tmp_I_E, _mm_sub_ps(_one, *_tmp_exp_E));
                
                
                _tmp_H_E  = *_tmp_E;
                _meanFR[i_node_vec] = _mm_add_ps(_meanFR[i_node_vec],_tmp_H_E);
                
                // Inhibitory population firing rate
                _tmp_I_I = _mm_sub_ps(_mm_mul_ps(_a_I,_mm_sub_ps(_mm_add_ps(_mm_add_ps(_w_I__I_0,_mm_mul_ps(_eIf_i, _ext_input[ext_inp_counter])),_mm_mul_ps(_J_NMDA, _S_i_E[i_node_vec])), _S_i_I[i_node_vec])),_b_I);
                ext_inp_counter++;
                *_tmp_exp_I   = _mm_mul_ps(_min_d_I, _tmp_I_I);
                tmp_exp_I[0]  = tmp_exp_I[0] != 0 ? expf(tmp_exp_I[0]) : 0.9;
                tmp_exp_I[1]  = tmp_exp_I[1] != 0 ? expf(tmp_exp_I[1]) : 0.9;
                tmp_exp_I[2]  = tmp_exp_I[2] != 0 ? expf(tmp_exp_I[2]) : 0.9;
                tmp_exp_I[3]  = tmp_exp_I[3] != 0 ? expf(tmp_exp_I[3]) : 0.9;
                *_tmp_I  = _mm_div_ps(_tmp_I_I, _mm_sub_ps(_one, *_tmp_exp_I));
                
                
                _tmp_H_I  = *_tmp_I;
                
                _r_i_E[i_node_vec] =  _tmp_H_E;
                _r_i_I[i_node_vec] =  _tmp_H_I;
                
                // Synaptic activity
                _S_i_I[i_node_vec] = _mm_add_ps(_S_i_I[i_node_vec],_mm_mul_ps(_dt,_mm_add_ps(_mm_mul_ps(_imintau_I, _S_i_I[i_node_vec]),_mm_mul_ps(_tmp_H_I,_gamma_I))));
                _S_i_E[i_node_vec] = _mm_add_ps(_S_i_E[i_node_vec],_mm_mul_ps(_dt, _mm_add_ps(_mm_mul_ps(_imintau_E, _S_i_E[i_node_vec]),_mm_mul_ps(_mm_mul_ps(_mm_sub_ps(_one, _S_i_E[i_node_vec]),_gamma),_tmp_H_E))));
            }
            for(j=0; j<nodes; j++){
                S_i_E[j] = S_i_E[j] < 1 ? S_i_E[j] : 1;
                S_i_E[j] = S_i_E[j] > 0 ? S_i_E[j] : 0;
                S_i_I[j] = S_i_I[j] < 1 ? S_i_I[j] : 1;
                S_i_I[j] = S_i_I[j] > 0 ? S_i_I[j] : 0;
            }
            
            ext_inp_counter -= nodes_vec;
            memcpy(&region_activity[ring_buf_pos], S_i_E, regions*sizeof( float ));
            // Shift ring-buff-pos
            ring_buf_pos = ring_buf_pos<(reg_act_size-regions) ? (ring_buf_pos+regions) : 0;
        }
        ext_inp_counter += nodes_vec;
        
        /*
         Compute BOLD for that time-step (subsampled to 1 ms)
         */
        for (j = 0; j < num_output_ts; j++) {
            // Excitatory populations
            bw_x_ex[j]  = bw_x_ex[j]  +  model_dt * (S_i_E[j] - kappa * bw_x_ex[j] - y * (bw_f_ex[j] - 1.0));
            f_tmp       = bw_f_ex[j]  +  model_dt * bw_x_ex[j];
            bw_nu_ex[j] = bw_nu_ex[j] +  model_dt * itau * (bw_f_ex[j] - powf(bw_nu_ex[j], ialpha));
            bw_q_ex[j]  = bw_q_ex[j]  +  model_dt * itau * (bw_f_ex[j] * (1.0 - powf(oneminrho,(1.0/bw_f_ex[j]))) / rho  - powf(bw_nu_ex[j],ialpha) * bw_q_ex[j] / bw_nu_ex[j]);
            bw_f_ex[j]  = f_tmp;
            
            /*
             // Inhibitory populations
             bw_x_in[j]  = bw_x_in[j]  +  model_dt * (S_i_I[j] - kappa * bw_x_in[j] - y * (bw_f_in[j] - 1.0));
             f_tmp       = bw_f_in[j]  +  model_dt * bw_x_in[j];
             bw_nu_in[j] = bw_nu_in[j] +  model_dt * itau * (bw_f_in[j] - powf(bw_nu_in[j], ialpha));
             bw_q_in[j]  = bw_q_in[j]  +  model_dt * itau * (bw_f_in[j] * (1.0 - powf(oneminrho,(1.0/bw_f_in[j]))) / rho  - powf(bw_nu_in[j],ialpha) * bw_q_in[j] / bw_nu_in[j]);
             bw_f_in[j]  = f_tmp;
             */
        }
        
        
        // Write out state variables as columns and time-steps as rows
        /*
        for (j = 0; j < num_output_ts; j++) {
            fprintf(Sout, "%.5f ",S_i_E[j]);
        }
        fprintf(Sout, "\n"); // End of output line
         */
        /*
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.4f ",S_i_I[j]);
        }
        
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.2f ",r_i_E[j]);
        }
        
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.2f ",r_i_I[j]);
        }
         */
        // divide by 10 to get sampling rate of 10 Hz
        for (i_node_vec = 0; i_node_vec < nodes_vec; i_node_vec++) {
            _EP_selfExc[i_node_vec]=_mm_mul_ps(_EP_selfExc[i_node_vec],  _pone);
            _EP_globInp[i_node_vec]=_mm_mul_ps(_EP_globInp[i_node_vec],  _pone);
            _EP_inhInp[i_node_vec] =_mm_mul_ps(_EP_inhInp[i_node_vec],  _pone);
            _EP_extInp[i_node_vec] =_mm_mul_ps(_EP_extInp[i_node_vec],  _pone);
            _EP_extInpinh[i_node_vec] =_mm_mul_ps(_EP_extInpinh[i_node_vec],  _pone);
        }
        /*
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.3f ",EP_selfExc[j]);
        }
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.3f ",EP_globInp[j]);
        }
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.3f ",EP_inhInp[j]);
        }
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.3f ",EP_extInp[j]);
        }
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.3f ",EP_extInpinh[j]);
        }
        */
        /*
        for (j = 0; j < num_output_ts; j++) {
            fprintf(WFout, "%.5f ",(EP_selfExc[j]+EP_globInp[j]+EP_extInp[j]));
        }
        fprintf(WFout, "\n"); // End of output line
         */
        // Excitatory populations' BOLD
//        for (j = 0; j < num_output_ts; j++) {
//            tmpBOLD = 100 / rho * V_0 * (k1 * (1 - bw_q_ex[j]) + k2 * (1 - bw_q_ex[j]/bw_nu_ex[j]) + k3 * (1 - bw_nu_ex[j]));
//            fprintf(WFout, "%.3f ",tmpBOLD);
//        }
        
        // Inhibitory populations' BOLD
        /*
         for (j = 0; j < num_output_ts; j++) {
         tmpBOLD = 100 / rho * V_0 * (k1 * (1 - bw_q_in[j]) + k2 * (1 - bw_q_in[j]/bw_nu_in[j]) + k3 * (1 - bw_nu_in[j]));
         fprintf(WFout, "%.3f ",tmpBOLD);
         }
         */

        
        
        if (ts_bold % TR_ts == 0) {
            BOLD_len_i++;
            
            for (j = 0; j < num_output_ts; j++) {
                BOLD_ex[j][BOLD_len_i] = 100 / rho * V_0 * (k1 * (1 - bw_q_ex[j]) + k2 * (1 - bw_q_ex[j]/bw_nu_ex[j]) + k3 * (1 - bw_nu_ex[j]));
                //BOLD_in[BOLD_len_i * num_output_ts + j] =  100 / rho * V_0 * (k1 * (1 - bw_q_ex[j]) + k2 * (1 - bw_q_ex[j]/bw_nu_ex[j]) + k3 * (1 - bw_nu_ex[j]));
            }
        }
        ts_bold++;
        
    }
    
    //fprintf(FCout, "%.10f %.10f %.10f %.10f %.10f %.2f \n\n", G,w_plus, J_NMDA, eIf_e, eIf_i, current_best_FRmean);
    for (i=0; i<(BOLD_len_i-BOLD_offset); i++) {
        for (j=0; j<num_output_ts; j++) {
            fprintf(FCout, "%.7f ",BOLD_ex[j][i+BOLD_offset]);
        }
        fprintf(FCout, "\n");
    }
    fprintf(FCout, "\n");
    fclose(FCout);
    //fclose(WFout);
    
//    _mm_free(n_conn_table);
//    _mm_free(region_activity);
//    _mm_free(n_conn_table_G_NMDA);
//    _mm_free(reg_globinp_p);
//    _mm_free(SC_cap);
//
    
    printf("%d finished. Execution took %.2f s\n", myid, (float)(time(NULL) - start));
    return 0;
}
