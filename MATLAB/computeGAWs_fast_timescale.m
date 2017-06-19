% Compute grand average waveforms (GAW) on fast time scale time-locked to 
% alpha phase fluctuation
%
% USAGE: 
% GAWnode = computeGAWnode_fast_timescale(waveform_file, node)
%
% INPUTS:
% waveform_file   - String that contains the filename of the ..txtwf
%                   waveform file output by the hybrid model C code
% node            - Scalar integer that specifies the brain network model
%                   node for which GAWnode should be computed
% OUTPUTS:
% GAWnode         - struct that contains grand average waveforms for node
%
% COMMENTS:
% This function sums and aggregates alpha cycles for one node. The actual
% grand average and grand standard deviations are then computed in the
% script computeGAWnode_fast_timescale_agg.m. This is necessary b/c: the data
% size is too large for a standard computer to have the GAW computed over
% all nodes in one run.


function GAWnode = computeGAWs_fast_timescale(waveform_file, node)

    % Read waveform file
    fid         =   fopen(waveform_file,'r');        
    data        =   fscanf(fid,'%f'); 
    fclose(fid);
    data        =   reshape(data,68*10,length(data)/(68*10))';


    %{
        Pre-allocate arrays
    %}
    num_guess               =   4500; % Initial guess
    
    % cycle length of alpha wave to get sharp GAW profiles
    cycle_lens              =   [95 105];
    
    % Arrays to store average
    GAWnode.grand_avg_alpha         = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_inhib         = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_FR_exc        = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_FR_inh        = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_Extern        = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_ExternInh     = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_Recurrent     = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_input         = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_inputInh      = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_ExcSyn        = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_InhSyn        = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_Global        = zeros(cycle_lens(2),1);
    GAWnode.grand_avg_BOLDexc       = zeros(cycle_lens(2),1);

    % Arrays to store standard deviation
    GAWnode.grand_avg_alphaSD       = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_inhibSD       = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_FR_excSD      = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_FR_inhSD      = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_ExternSD      = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_ExternInhSD   = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_RecurrentSD   = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_inputSD       = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_inputInhSD    = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_ExcSynSD      = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_InhSynSD      = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_GlobalSD      = zeros(num_guess,cycle_lens(2));
    GAWnode.grand_avg_BOLDexcSD     = zeros(num_guess,cycle_lens(2));

    % Binned alpha power histogramm
    GAWnode.bin1                    = zeros(num_guess,1);
    GAWnode.bin2                    = zeros(num_guess,1);
    GAWnode.bin3                    = zeros(num_guess,1);
    GAWnode.bin4                    = zeros(num_guess,1);
    GAWnode.bin5                    = zeros(num_guess,1);
    GAWnode.bin6                    = zeros(num_guess,1);
    
    
    %{
        Extract full waveforms
    %}

    % Butterworth band-pass filter in alpha range
    [b,a]                   =   butter(1, [8 12]/(1000/2));
    
    % Extract and filter injected input
    extern_inp              =   data(:,node+68*7);
    extern_inp_filt         =   filtfilt(b,a,extern_inp);
    
    % Compute instantaneous amplitude (ongoing alpha power)
    extern_inp_filt_hilbert =   hilbert(extern_inp_filt);
    inst_ampl               =   abs(extern_inp_filt_hilbert);
    GAWnode.alphapower      =   inst_ampl;

    % Firing rate of exc. pop.
    FR_exc                  =   data(:,node+68*2);
    % Firing rate of inh. pop.
    FR_inh                  =   data(:,node+68*3);
    % Synaptic gating of exc. pop.
    Syn_exc                 =   data(:,node);
    % Synaptic gating of inh. pop.
    Syn_inh                 =   data(:,node+68);
    % Injected input exc. pop.
    externinp               =   data(:,node+68*7);
    % Injected input inh. pop. (scaled version of inj. inp. to exc. pop.)
    externinpinh            =   data(:,node+68*8);
    % Long-range input
    Global_inp              =   data(:,node+68*5);
    % Feedback inhibition
    Inhib_inp               =   data(:,node+68*6);
    % Recurrent excitation
    Exc_inp                 =   data(:,node+68*4);
    % Sum of all inputs to exc. pop.
    Input_ExcPop            =   externinp + 0.382 - Inhib_inp + Global_inp;
    % Sum of all inputs to inh. pop.
    Input_InhPop            =   0.7*0.382 - Syn_inh + externinpinh;
    % fMRI
    bold_exc                =   data(:,node+68*9);

 
    
    %{
        Compute grand average waveforms time-locked to alpha phase fluctuation
    %}

    counter                 =   0;
    previous_cycle_start    =   -1;
    
    % Scan through alpha wave 
    % Edges are ignored to avoid edge artifacts
    for ii = 1000*50*1.94:length(extern_inp_filt)-1000*50*1.94,
        
        % Zero-crossing of alpha wave marks the beginning of an alpha wave
        if extern_inp_filt(ii-1) < 0 && extern_inp_filt(ii) >= 0, 
            if previous_cycle_start == -1,
                previous_cycle_start = ii;            
            else        
                % Compute length of current alpha cycle
                cycle_length = ii-1 - previous_cycle_start + 1;

                % Only use for GAW computation if length is in prescribed
                % boundaries (to get sharp GAW profiles)
                if cycle_length > cycle_lens(1) && cycle_length < cycle_lens(2)                

                    counter=counter+1;

                    % Mean firing rate of current cycle
                    avg_FR_curr_cycle   = mean(FR_exc(previous_cycle_start:(ii-1)));  
                    
                    % Divide cycle into six sub-epochs, compute mean firing
                    % rate during each epoch
                    GAWnode.bin1(counter)       = mean(FR_exc(previous_cycle_start: (previous_cycle_start + round(cycle_length/6) - 1)))./avg_FR_curr_cycle;
                    GAWnode.bin2(counter)       = mean(FR_exc((previous_cycle_start + round(1*cycle_length/6)) : (previous_cycle_start + round(2*cycle_length/6) - 1)))./avg_FR_curr_cycle;
                    GAWnode.bin3(counter)       = mean(FR_exc((previous_cycle_start + round(2*cycle_length/6)) : (previous_cycle_start + round(3*cycle_length/6) - 1)))./avg_FR_curr_cycle;
                    GAWnode.bin4(counter)       = mean(FR_exc((previous_cycle_start + round(3*cycle_length/6)) : (previous_cycle_start + round(4*cycle_length/6) - 1)))./avg_FR_curr_cycle;
                    GAWnode.bin5(counter)       = mean(FR_exc((previous_cycle_start + round(4*cycle_length/6)) : (previous_cycle_start + round(5*cycle_length/6) - 1)))./avg_FR_curr_cycle;
                    GAWnode.bin6(counter)       = mean(FR_exc((previous_cycle_start + round(5*cycle_length/6)) : (previous_cycle_start + round(6*cycle_length/6) - 1)))./avg_FR_curr_cycle;

                    % Sum waveforms to compute average for each time step
                    % later (during aggregation step)
                    GAWnode.grand_avg_alpha     = GAWnode.grand_avg_alpha       + extern_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_inhib     = GAWnode.grand_avg_inhib       - Inhib_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_InhSyn    = GAWnode.grand_avg_InhSyn      + Syn_inh(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_FR_exc    = GAWnode.grand_avg_FR_exc      + FR_exc(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_FR_inh    = GAWnode.grand_avg_FR_inh      + FR_inh(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_Extern    = GAWnode.grand_avg_Extern      + externinp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_ExternInh = GAWnode.grand_avg_ExternInh   + externinpinh(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_Recurrent = GAWnode.grand_avg_Recurrent   + Exc_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);                
                    GAWnode.grand_avg_input     = GAWnode.grand_avg_input       + Input_ExcPop(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);        
                    GAWnode.grand_avg_inputInh  = GAWnode.grand_avg_inputInh    + Input_InhPop(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);                   
                    GAWnode.grand_avg_ExcSyn    = GAWnode.grand_avg_ExcSyn      + Syn_exc(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_Global    = GAWnode.grand_avg_Global      + Global_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_BOLDexc   = GAWnode.grand_avg_BOLDexc     + bold_exc(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);

                    % Store waveforms to compute SD for each time step
                    % later (during aggregation step)
                    GAWnode.grand_avg_alphaSD(counter,:)        = extern_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_inhibSD(counter,:)        = -Inhib_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_FR_excSD(counter,:)       = FR_exc(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_FR_inhSD(counter,:)       = FR_inh(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_ExternSD(counter,:)       = externinp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_ExternInhSD(counter,:)    = externinpinh(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_RecurrentSD(counter,:)    = Exc_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_inputSD(counter,:)        = Input_ExcPop(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1); 
                    GAWnode.grand_avg_inputInhSD(counter,:)     = Input_InhPop(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);            
                    GAWnode.grand_avg_ExcSynSD(counter,:)       = Syn_exc(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_InhSynSD(counter,:)       = Syn_inh(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_GlobalSD(counter,:)       = Global_inp(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                    GAWnode.grand_avg_BOLDexcSD(counter,:)      = bold_exc(previous_cycle_start:previous_cycle_start+cycle_lens(2)-1);
                end

                % Time step of the end/beginning of this/next cycle 
                previous_cycle_start = ii; 
            end
        end 
    end

    % Store results in output struct (cut off empty parts of arrays)
    GAWnode.bin1                       =   bin1(1:counter);
    GAWnode.bin2                       =   bin2(1:counter);
    GAWnode.bin3                       =   bin3(1:counter);
    GAWnode.bin4                       =   bin4(1:counter);
    GAWnode.bin5                       =   bin5(1:counter);
    GAWnode.bin6                       =   bin6(1:counter);
    GAWnode.grand_avg_alphaSD          =   GAWnode.grand_avg_alphaSD(1:counter,:);
    GAWnode.grand_avg_inhibSD          =   GAWnode.grand_avg_inhibSD(1:counter,:);
    GAWnode.grand_avg_FRSD             =   GAWnode.grand_avg_FR_excSD(1:counter,:);
    GAWnode.grand_avg_ExternSD         =   GAWnode.grand_avg_ExternSD(1:counter,:);
    GAWnode.grand_avg_ExternInhSD      =   GAWnode.grand_avg_ExternInhSD(1:counter,:);
    GAWnode.grand_avg_RecurrentSD      =   GAWnode.grand_avg_RecurrentSD(1:counter,:);
    GAWnode.grand_avg_inputSD          =   GAWnode.grand_avg_inputSD(1:counter,:);
    GAWnode.grand_avg_inputInhSD       =   GAWnode.grand_avg_inputInhSD(1:counter,:);
    GAWnode.grand_avg_ExcSynSD         =   GAWnode.grand_avg_ExcSynSD(1:counter,:);
    GAWnode.grand_avg_GlobalSD         =   GAWnode.grand_avg_GlobalSD(1:counter,:);
    GAWnode.grand_avg_InhSynSD         =   GAWnode.grand_avg_InhSynSD(1:counter,:);
    GAWnode.grand_avg_FR_inhSD         =   GAWnode.grand_avg_FR_inhSD(1:counter,:);
    GAWnode.grand_avg_BOLDexcSD        =   GAWnode.grand_avg_BOLDexcSD(1:counter,:);
    GAWnode.counter                    =   counter;

end


