% Compute alpha regressor from EEG source activity
%
% USAGE: 
% alpha_reg = compute_alpha_regressor(source_activity, regionsMap, HRF_resEEG)
%
% INPUTS:
% source_activity - [68 x 259184] matrix containing EEG source activity for 
%                   68 regions and 259184 time points (200 Hz sampling rate)
%                   (data/empirical_source_activity/source_activity.mat)
% regionsMap      - [68 x 1] vector that contains the region sorting of source
%                   activity as outputted by Brainstorm 
%                   (data/regionsMap.mat)
% HRF_resEEG      - [1 x 6498] vector that contains the hemodynamic
%                   response function sampled at 200 Hz
%                   (data/HRF_200Hz.mat)
%
% OUTPUTS:
% alpha_reg       - struct that contains alpha regressor and filtered
%                   alpha regressor (edges are discarded due to edge 
%                   effects from filtering and from convolution with HRF)



function alpha_reg = compute_alpha_regressor(source_activity, regionsMap, HRF_resEEG)
    % Sorting of Desikan-Killiany atlas regions in SC matrices
    SCmat_sorting=[1001:1003,1005:1035,2001:2003,2005:2035];
    
    % Generate butterworth filter for alpha range and for resting-state
    % slow oscillations range
    [b_hi,a_hi]     =   butter(1, [8 12]/(200/2));
    [b_lo,a_lo]     =   butter(1, [0.1]/((1/1.94)/2));


    % Initialize output arrays (shorter than full fMRI time series to
    % discard edge effects from filtering and convolution with HRF)
    alpha_reg       =   zeros(651,68);    
    alpha_reg_filt  =   zeros(651,68); 
    
    % Iterate over regions
    for ii = 1:68,
        % Get SC matrix sorting
        regindSAC           =   find(regionsMap==SCmat_sorting(ii));
        region_ts           =   source_activity(regindSAC,:);
        
        % Filter in alpha range
        region_ts_filt      =   filtfilt(b_hi,a_hi,region_ts);
        
        % Hilbert transform to get instantaneous amplitude
        region_ts_filt_hilb =   hilbert(region_ts_filt);
        inst_ampl           =   abs(region_ts_filt_hilb);
        
        % Convolution with HRF
        inst_ampl_HRF       =   conv(inst_ampl(100:end-100),HRF_resEEG,'valid');
        
        % Downsample to BOLD sampling rate (TR = 1.94 s)
        alpha_reg(:,ii)     =   downsample(inst_ampl_HRF(100:end-100),388);
        alpha_reg_filt(:,ii)=   filtfilt(b_lo,a_lo,alpha_reg(:,ii));
    end

    % Fill output struct
    alpha_reg.alpha_reg         = alpha_reg;
    alpha_reg.alpha_reg_filt    = alpha_reg_filt;

end
