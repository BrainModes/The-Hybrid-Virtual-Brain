% Generate hybrid model input from EEG source activity
%
% USAGE: 
% preprocess_external_input(output_fileid, source_activity, regionsMap)
%
% INPUTS:
% output_fileid   - String that contains the filestem for the output textfile
% source_activity - [68 x 259184] matrix containing EEG source activity for 
%                   68 regions and 259184 time points (200 Hz sampling rate)
%                   (stored in data/empirical_source_activity/source_activity.mat)
% regionsMap      - [68 x 1] vector that contains the region sorting of source
%                   activity as outputted by Brainstorm 
%
%
% OUTPUTS:
% [output_fileid '.txt']            - hybrid model input activity
% [output_fileid '_randperm.txt']   - random permutation of input activity
%
% COMMENTS:
% Region-wise source imaging results from Brainstorm are re-sorted to match
% the region-sorting of SC matrices and upsampled to 1000 Hz

function generate_hybrid_model_input(output_fileid, source_activity, regionsMap)

    % Sorting of Desikan-Killiany atlas regions in SC matrices
    SCmat_sorting   =   [1001:1003,1005:1035,2001:2003,2005:2035];

    % Sampling rate of model input activity (Hz)
    srate           =   1000;

    % Number of TRs in source activity (1 TR = 1.94 s)
    numscans        =   668;

    % Original and interpolated time points
    x               =   1:(srate/200):(srate*1.94*numscans);
    xi              =   1:(srate*1.94*numscans);

    % Initialize output arrays
    sorted_sac      =   zeros(length(xi),68);
    sorted_sacrand  =   zeros(length(xi),68);

    % Re-sort and upsample source activity; generate random permutation
    for ii = 1:68,
        regindSAC           =   find(regionsMap==SCmat_sorting(ii));
        ts                  =   zscore(source_activity(regindSAC,:));
        sorted_sac(:,ii)    =   interp1(x,ts,xi,'spline','extrap');
        sorted_sacrand(:,ii)=   sorted_sac(randperm(size(sorted_sac,1)),ii);
    end

    % Write output
    dlmwrite([output_fileid '.txt'],sorted_sac,'delimiter',' ','precision','%.4f');
    dlmwrite([output_fileid '_randperm.txt'],sorted_sac,'delimiter',' ','precision','%.4f');

end
