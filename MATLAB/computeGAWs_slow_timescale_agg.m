% Aggregates grand average waveform (GAW) results over all nodes and 
% subjects for activity on a slow time scale (< 0.1 Hz)
%
% USAGE: 
% GAW = computeGAWs_fast_timescale_agg(GAWnode_filenames)
%
% INPUTS:
% GAWnode_filenames - cell array that contains filenames for results from
%                     extract_statevars.m (68 regions for 15 subjects)
%
% OUTPUTS:
% GAW               - grand average waveform over all nodes and subjects



function GAW = computeGAWs_slow_timescale_agg(GAWnode_filenames)

    % Initialize arrays
    GAW.EPSC    = zeros(129398,1);
    GAW.EPSCinh = zeros(129398,1);
    GAW.SynINH  = zeros(129398,1);
    GAW.SynEXC  = zeros(129398,1);
    GAW.IPSC    = zeros(129398,1);
    GAW.FRexc   = zeros(129398,1);
    GAW.FRinh   = zeros(129398,1);
    GAW.GlobaI  = zeros(129398,1);
    GAW.BOLD    = zeros(129398,1);
    GAW.INPexc  = zeros(129398,1);
    GAW.INPinh  = zeros(129398,1);

    % Iterate over all GAW result files and compute sum of waveforms
    for ii = length(GAWnode_filenames)
        load(GAWnode_filenames{ii});
        GAW.EPSC    = GAW.EPSC    + sum(statevars.externinp')';
        GAW.EPSCinh = GAW.EPSCinh + sum(statevars.externinpinh')';
        GAW.SynINH  = GAW.SynINH  + sum(statevars.Syn_inh_ts')';
        GAW.SynEXC  = GAW.SynEXC  + sum(statevars.Syn_exc')';
        GAW.IPSC    = GAW.IPSC    + sum(statevars.Inhib_inp')';
        GAW.FRexc   = GAW.FRexc   + sum(statevars.FR_exc')';
        GAW.FRinh   = GAW.FRinh   + sum(statevars.FR_inh')';
        GAW.GlobaI  = GAW.GlobaI  + sum(statevars.Global_inp')';
        GAW.BOLD    = GAW.BOLD    + sum(statevars.bold_exc')';
        GAW.INPexc  = GAW.INPexc  + sum(statevars.Input_ExcPop')';
        GAW.INPinh  = GAW.INPinh  + sum(statevars.Input_InhPop')';
    end

    % Compute average waveform
    GAW.EPSC    = GAW.EPSC    / (68*15) ;
    GAW.EPSCinh = GAW.EPSCinh / (68*15) ;
    GAW.SynINH  = GAW.SynINH  / (68*15) ;
    GAW.SynEXC  = GAW.SynEXC  / (68*15) ;
    GAW.IPSC    = GAW.-1*IPSC / (68*15) ;
    GAW.FRexc   = GAW.FRexc   / (68*15) ;
    GAW.FRinh   = GAW.FRinh   / (68*15) ;
    GAW.GlobaI  = GAW.GlobaI  / (68*15) ;
    GAW.BOLD    = GAW.BOLD    / (68*15) ;
    GAW.INPexc  = GAW.INPexc  / (68*15) ;
    GAW.INPinh  = GAW.INPinh  / (68*15) ;

    % Compute moving averages
    movavg_window       =   50;
    GAW.movavgFR_exc    =   nan(size(FRexc,1),1);
    GAW.movavgFR_inh    =   nan(size(FRinh,1),1);
    for p = (movavg_window + 1) : size(FRexc,1) - movavg_window,
        GAW.movavgFR_exc(p) = mean(FRexc(p-movavg_window:p+movavg_window));
        GAW.movavgFR_inh(p) = mean(FRinh(p-movavg_window:p+movavg_window));
    end

end
