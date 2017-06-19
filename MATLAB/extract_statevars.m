% Extract and downsample state variables from simulation results
%
% USAGE: 
% statevars = extract_statevars(waveform_file)
%
% INPUTS:
% waveform_file   - String that contains the filename of the ..txtwf
%                   waveform file output by the hybrid model C code
%
% OUTPUTS:
% statevars       - struct that contains state variables of the hybrid 
%                   model



function statevars = extract_statevars(waveform_file)

    % Read waveform file
    fid         =   fopen(waveform_file,'r');        
    data        =   fscanf(fid,'%f'); 
    fclose(fid);
    data        =   reshape(data,68*10,length(data)/(68*10))';

    %{
    Downsample to 100 Hz, fill output structure
    %}
    nodes                     = (1:68);
    
    % Firing rate of excitatory population
    statevars.FR_exc        = downsample(data(:,nodes+68*2),10);
    
    % Firing rate of inhibitory population
    statevars.FR_inh        = downsample(data(:,nodes+68*3),10);
    
    % Synaptic gating of excitatory population
    statevars.Syn_exc       = downsample(data(:,nodes),10);
    
    % Synaptic gating of inhibitory population
    statevars.Syn_inh       = downsample(data(:,nodes+68),10);
    
    % Source activity injected into exc. pop.
    statevars.externinp     = downsample(data(:,nodes+68*7),10);
    
    % Source activity injected into inh. pop. (different scaling than into
    % exc. pop.)
    statevars.externinpinh  = downsample(data(:,nodes+68*8),10);
    
    % Long-range input 
    statevars.Global_inp    = downsample(data(:,nodes+68*5),10);
    
    % Feedback inhibition
    statevars.Inhib_inp     = downsample(data(:,nodes+68*6),10);
    
    % Recurrent excitation
    statevars.Exc_inp       = downsample(data(:,nodes+68*4),10);
    
    % Sum of all inputs to exc. pop.
    statevars.Input_ExcPop  = statevars.unfexterninp + 0.382 - statevars.unfInhib_inp + statevars.unfGlobal_inp;
    
    % Sum of all inputs to inh. pop.
    statevars.Input_InhPop  = 0.7*0.382 - statevars.unfSyn_inh_ts + statevars.unfexterninpinh;
    
    % fMRI of exc. pop.
    statevars.bold_exc      = downsample(data(:,nodes+68*9),10);

end


