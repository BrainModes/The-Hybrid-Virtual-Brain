# The-Hybrid-Virtual-Brain
Quick instructions

1. Compile the C code with a compiler with support for SSE instructions enabled. I found the following combination of flags yields good performance (in terms of simulation speed) with the GNU C compiler

gcc  -Wall -std=c99 -msse2 -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -m64 -lm TVBii.c -o tvb


2. Create a folder “input” and a folder “output” in the folder where the newly created program binary is stored.


3. Create input files for each subject-specific brain model using the MATLAB script 'generate_input_SC.m' and store them in the folder 'input'. For each subject three input files are required. They contain (i) connection strengths, (ii) connection lengths (unit: mm), and (iii) the node-ids of the connected nodes in files (i) and (ii). The MATLAB script converts structural connectivity from a matrix representation (strengths matrix and delay matrix) into a list format that can be read by the hybrid model program. In addition, the user needs to supply files that contain the time series (sampled at 1000 Hz) that are to be injected into the model populations. These files must also be stored in the folder 'input'. Different types of input were considered during the experimental phase of the project (EEG source activity, randomly permutated EEG source activity inputs, temporal differences of source activity, artificial alpha, etc.), for an overwiew see the different cases beginning at line 239 of TVBii.c. The different modes can be switched with the <input_case> argument shown in the next step. The standard case for source activity injection is case "0" and it expects a file of the format <subject\_id>\_fullextinp1000.txt in the 'input' folder.


4. Run the program in the folder from a terminal. The program expects three arguments:

./tvb <parameter_file> <subject_id> <input_case>

<parameter_file>: filename of a simple text file (stored in the same folder as <program>) that contains a sequence of three floating point numbers that specify (in this order)
* long-range coupling strength parameter G
* inhibitory population input scaling parameter w_BG^I
* and the ratio of inhibitory to excitatory input strengths w_BG^I / w_BG^E  (“optimal” ratios w_BG^I / w_BG^E found during parameter space exploration were >1)

<subject_id>: a string that specifies the prefix of all subject-specific input files (stored in folder ‘./input’ relative to program folder)

<input_case>: switches between eight different input files and output folders, see lines 239f of the code. It basically just adds a suffix to the input filename output path folder, for example

'input/<subject_name>_fullextinp1000.txt'


5. Outputs are stored in the folder 'output'. Files with prefix 'BOLD_' contain simulated fMRI time series (time x nodes)- Files with prefix 'SV_' contain state variables (time x n * state_variable). Please refer to lines 885 to 923 of the C code for the sorting of state variables. If you perform a large number of simulations, consider commenting out these lines, as the produced state variable files will be very large (several GB per file depending on simulated time).


Questions? michael.schirner@charite.de or petra.ritter@charite.de
