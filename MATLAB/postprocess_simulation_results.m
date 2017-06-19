% Reads simulation results, computes correlation between empirical and
% simulated time series and functional connectivities
%
% USAGE: 
% sim_res = postprocess_simulation_results(simulated_fMRI, empirical_fMRI)
%
% INPUTS:
% simulated_fMRI    - string that contains filename of hybrid model
%                     fMRI simulation results for subject
% empirical_fMRI    - [N,M] array containing N time points and M regions of
%                     empirical fMRI activity for subject
%                     (empirical_fMRI.mat)
% RSN               - [N,K] array containing N time points and K empirical 
%                     RSN temporal modes for subject (GroupICA_results.mat)
% RSNroi            - struct that contains mapping between RSNs and BNM
%                     regions (RSN_to_region.mat)
%
% OUTPUTS:
% sim_res           - struct that contains correlations between simulated
%                     and empirical time series and functional connectivity

function sim_res = postprocess_simulation_results(simulated_fMRI, empirical_fMRI, RSN, RSNroi)

    % logical array to map RSNs to BNM regions
    RSNidx_loose = zeros(68,10);
    for ii = 1:10
        RSNidx_loose(RSNroi(ii).SCrois_loose,ii) = 1;    
    end
    RSNidx_loose = logical(RSNidx_loose);


    % Compute low-pass filtered versions of fMRI/RSN activity
    [b,a]       =       butter(1, [0.1]/((1/1.94)/2));
    for ii = 1:68,
        empirical_fMRI_filt(:,ii)     =   filtfilt(b,a,empirical_fMRI(:,ii));
    end
    for ii = 1:10,
        RSN_filt(:,ii)      =   filtfilt(b,a,RSN(:,ii));
    end

    % Compute functional connectivity
    static_FCemp         = corr(empirical_fMRI(1+5+7:640+5+7,:));
    static_FCemp_filt    = corr(empirical_fMRI_filt(1+5+7:640+5+7,:));
    Isubdiag = find(tril(ones(68),-1)); % Indexes of all the values below the diagonal.


    % Read simulated fMRI
    fid                     =   fopen(simulated_fMRI,'r');
    data                    =   fscanf(fid,'%f');
    fclose(fid);


    % Pre-allocate result arrays
    ts_cc                   =   zeros(length(data)/(6+655*68),12);
    ts_cc_filt              =   zeros(length(data)/(6+655*68),12);
    static_FC_cc            =   zeros(length(data)/(6+655*68),2);
    static_FC_cc_filt       =   zeros(length(data)/(6+655*68),2);


    % Loop over all parameter sets
    for ii = 0:length(data)/(6+655*68)-1,   

        % Store simulation parameters / metadata
        ts_cc(ii+1,1:6)     =   data(1+ii*6+ii*655*68:6+ii*6+ii*655*68);

        % Extract simulated fMRI for parameter set
        sim_ts              =   reshape(data(7+ii*6+ii*655*68:7+ii*6+(ii+1)*655*68-1),68,655)';

        % Compute low-pass filtered version of simulated fMRI
        for node = 1:68,
            simts_filt(:,node)    =   filtfilt(b,a,sim_ts(:,node));
        end

        % Compute average fMRI time series correlations over all 68 simulated and 
        % empirical region-pairs for time shifts to account for lag of HRF.
        % Leave offset for edge effects of filtering.    
        clear tmp_res
        for shift = 0:15
            tmp_tmp_res         =   diag(corr(sim_ts(6:645,:),empirical_fMRI(6+shift:645+shift,:)));
            tmp_res(shift+1)    =   mean(tmp_tmp_res(~isnan(tmp_tmp_res)));
        end
        [ts_cc(ii+1,7), ts_cc(ii+1,8)]    =   max(tmp_res);

        % same for filtered data
        clear tmp_res
        for shift = 0:15
            tmp_tmp_res         =   diag(corr(simts_filt(6:645,:),empirical_fMRI_filt(6+shift:645+shift,:)));
            tmp_res(shift+1)    =   mean(tmp_tmp_res(~isnan(tmp_tmp_res)));
        end
        [ts_cc_filt(ii+1,7), ts_cc_filt(ii+1,8)]    =   max(tmp_res);


        % Compute average RSN time series correlations using all correlations 
        % between all 9 empirical RSNs and their corresponding simulated
        % regions (cerebellar RSN is excluded)
        clear tmp_res
        for shift = 0:15
            tmp_tmp_res         =   corr(sim_ts(1+5:640+5,:),RSN(1+shift+5:640+shift+5,:));
            tmp_res(shift+1,:)  =   [mean(tmp_tmp_res(RSNidx_loose(:,1),1)) mean(tmp_tmp_res(RSNidx_loose(:,2),2)) mean(tmp_tmp_res(RSNidx_loose(:,3),3)) mean(tmp_tmp_res(RSNidx_loose(:,4),4)) mean(tmp_tmp_res(RSNidx_loose(:,6),6)) mean(tmp_tmp_res(RSNidx_loose(:,7),7)) mean(tmp_tmp_res(RSNidx_loose(:,8),8)) mean(tmp_tmp_res(RSNidx_loose(:,9),9)) mean(tmp_tmp_res(RSNidx_loose(:,10),10))];
        end
        ts_cc(ii+1,9)           =   max(tmp_res(:));
        ts_cc(ii+1,10)          =   max(mean(tmp_res'));      

        % same for filtered data
        clear tmp_res
        for shift = 0:15
            tmp_tmp_res         =   corr(simts_filt(1+5:640+5,:),RSN_filt(1+shift+5:640+shift+5,:));
            tmp_res(shift+1,:)  =   [mean(tmp_tmp_res(RSNidx_loose(:,1),1)) mean(tmp_tmp_res(RSNidx_loose(:,2),2)) mean(tmp_tmp_res(RSNidx_loose(:,3),3)) mean(tmp_tmp_res(RSNidx_loose(:,4),4)) mean(tmp_tmp_res(RSNidx_loose(:,6),6)) mean(tmp_tmp_res(RSNidx_loose(:,7),7)) mean(tmp_tmp_res(RSNidx_loose(:,8),8)) mean(tmp_tmp_res(RSNidx_loose(:,9),9)) mean(tmp_tmp_res(RSNidx_loose(:,10),10))];
        end
        ts_cc_filt(ii+1,9)      =   max(tmp_res(:));
        ts_cc_filt(ii+1,10)     =   max(mean(tmp_res')); 


        % Correlation between simulated and empirical FC
        static_FC_sim                       =   corr(sim_ts(1+5:640+5,:));        
        static_FC_cc(ii+1,1)                =   corr(static_FC_sim(Isubdiag),static_FCemp(Isubdiag));

        % same for filtered data
        static_FC_sim_filt                  =   corr(simts_filt(1+5:640+5,:));        
        static_FC_cc_filt(ii+1,1)           =   corr(static_FC_sim_filt(Isubdiag),static_FCemp_filt(Isubdiag));

    end

    sim_res.ts_cc               = ts_cc;
    sim_res.ts_cc_filt          = ts_cc_filt;
    sim_res.static_FC_cc        = static_FC_cc;
    sim_res.static_FC_cc_filt   = static_FC_cc_filt;


end
