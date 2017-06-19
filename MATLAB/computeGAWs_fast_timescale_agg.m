% Aggregates grand average waveform (GAW) results over all nodes and 
% subjects
%
% USAGE: 
% GAW = computeGAWs_fast_timescale_agg(GAWnode_filenames)
%
% INPUTS:
% GAWnode_filenames - cell array that contains filenames for results from
%                     computeGAWs_fast_timescale.m. 
%
% OUTPUTS:
% GAW               - grand average waveform over all nodes and subjects
%
% COMMENTS:
% This function computes the average over all individual GAWnode results
% for each subject and node output by computeGAWs_fast_timescale.m.

function GAW = computeGAWs_fast_timescale_agg(GAWnode_filenames)

    % Initialize arrays
    GAW.grand_avg_alphaSDc          = zeros(1,105); 
    GAW.grand_avg_inhibSDc          = zeros(1,105); 
    GAW.grand_avg_FRSDc             = zeros(1,105); 
    GAW.grand_avg_ExternSDc         = zeros(1,105); 
    GAW.grand_avg_ExternInhSDc		= zeros(1,105); 
    GAW.grand_avg_RecurrentSDc      = zeros(1,105); 
    GAW.grand_avg_inputSDc  		= zeros(1,105); 
    GAW.grand_avg_inputInhSDc		= zeros(1,105); 
    GAW.grand_avg_ExcSynSDc 		= zeros(1,105); 
    GAW.grand_avg_GlobalSDc 		= zeros(1,105); 
    GAW.grand_avg_InhSynSDc 		= zeros(1,105); 
    GAW.grand_avg_FR_inhSDc 		= zeros(1,105); 
    GAW.grand_avg_BOLDexcSDc		= zeros(1,105); 
    GAW.bin1                        = zeros(1,105); 
    GAW.bin2                        = zeros(1,105); 
    GAW.bin3                        = zeros(1,105); 
    GAW.bin4                        = zeros(1,105); 
    GAW.bin5                        = zeros(1,105); 
    GAW.bin6                        = zeros(1,105); 
    SD_GAW.grand_avg_alphaSDc		= zeros(1,105); 
    SD_GAW.grand_avg_inhibSDc 		= zeros(1,105); 
    SD_GAW.grand_avg_FRSDc			= zeros(1,105); 
    SD_GAW.grand_avg_ExternSDc		= zeros(1,105); 
    SD_GAW.grand_avg_ExternInhSDc	= zeros(1,105); 
    SD_GAW.grand_avg_RecurrentSDc 	= zeros(1,105); 
    SD_GAW.grand_avg_inputSDc  		= zeros(1,105); 
    SD_GAW.grand_avg_inputInhSDc	= zeros(1,105); 
    SD_GAW.grand_avg_ExcSynSDc 		= zeros(1,105); 
    SD_GAW.grand_avg_GlobalSDc 		= zeros(1,105); 
    SD_GAW.grand_avg_InhSynSDc 		= zeros(1,105); 
    SD_GAW.grand_avg_FR_inhSDc 		= zeros(1,105); 
    SD_GAW.grand_avg_BOLDexcSDc		= zeros(1,105); 
    SD_GAW.bin1                     = zeros(1,105); 
    SD_GAW.bin2                     = zeros(1,105); 
    SD_GAW.bin3                     = zeros(1,105); 
    SD_GAW.bin4                     = zeros(1,105); 
    SD_GAW.bin5                     = zeros(1,105); 
    SD_GAW.bin6                     = zeros(1,105); 

    % Initialize element counters for each bin
    lenb1 = 0;
    lenb2 = 0;
    lenb3 = 0;
    lenb4 = 0;
    lenb5 = 0;    
    lenb6 = 0;

    % Sum up all node-wise GAW results
    total_alpha_cycles = 0;
    for ii = 1:length(GAWnode_filenames)

        % Load GAWnode file
        load(GAWnode_filenames{ii});

        % Sum waveforms to compute average waveform    
        GAW.grand_avg_alphaSDc          = GAW.grand_avg_alphaSDc        +   sum(    GAWnode.grand_avg_alphaSDc    );
        GAW.grand_avg_inhibSDc          = GAW.grand_avg_inhibSDc        +   sum(    GAWnode.grand_avg_inhibSDc    );
        GAW.grand_avg_FRSDc             = GAW.grand_avg_FRSDc           +   sum(    GAWnode.grand_avg_FRSDc       );
        GAW.grand_avg_ExternSDc         = GAW.grand_avg_ExternSDc       +   sum(    GAWnode.grand_avg_ExternSDc   );
        GAW.grand_avg_ExternInhSDc		= GAW.grand_avg_ExternInhSDc	+   sum(    GAWnode.grand_avg_ExternInhSDc);
        GAW.grand_avg_RecurrentSDc      = GAW.grand_avg_RecurrentSDc	+   sum(    GAWnode.grand_avg_RecurrentSDc);
        GAW.grand_avg_inputSDc  		= GAW.grand_avg_inputSDc        +   sum(    GAWnode.grand_avg_inputSDc    );
        GAW.grand_avg_inputInhSDc		= GAW.grand_avg_inputInhSDc     +   sum(    GAWnode.grand_avg_inputInhSDc );
        GAW.grand_avg_ExcSynSDc 		= GAW.grand_avg_ExcSynSDc       +   sum(    GAWnode.grand_avg_ExcSynSDc   );
        GAW.grand_avg_GlobalSDc 		= GAW.grand_avg_GlobalSDc       +   sum(    GAWnode.grand_avg_GlobalSDc   );
        GAW.grand_avg_InhSynSDc 		= GAW.grand_avg_InhSynSDc       +   sum(    GAWnode.grand_avg_InhSynSDc   );
        GAW.grand_avg_FR_inhSDc 		= GAW.grand_avg_FR_inhSDc       +   sum(    GAWnode.grand_avg_FR_inhSDc   );
        GAW.grand_avg_BOLDexcSDc		= GAW.grand_avg_BOLDexcSDc      +   sum(    GAWnode.grand_avg_BOLDexcSDc  );

        % Sum up number of elements for each bin
        lenb1                           = lenb1     +   length( GAWnode.bin1   );
        lenb2                           = lenb2     +   length( GAWnode.bin2   );
        lenb3                           = lenb3     +   length( GAWnode.bin3   );
        lenb4                           = lenb4     +   length( GAWnode.bin4   );
        lenb5                           = lenb5     +   length( GAWnode.bin5   );    
        lenb6                           = lenb6     +   length( GAWnode.bin6   );

        % Sum up firing rates in each bin
        GAW.bin1                        = GAW.bin1  +   sum(    GAWnode.bin1   );
        GAW.bin2                        = GAW.bin2  +   sum(    GAWnode.bin2   );
        GAW.bin3                        = GAW.bin3	+   sum(    GAWnode.bin3   );
        GAW.bin4                        = GAW.bin4	+   sum(    GAWnode.bin4   );
        GAW.bin5                        = GAW.bin5	+   sum(    GAWnode.bin5   );
        GAW.bin6                        = GAW.bin6	+   sum(    GAWnode.bin6   );

        % Sum up number of alpha cycles
        total_alpha_cycles              = total_alpha_cycles + GAWnode.counter ;
    end


    % Compute average waveforms by dividing each time step by number of alpha
    % cycles
    GAW.grand_avg_alphaSDc          =  GAW.grand_avg_alphaSDc       / total_alpha_cycles; 
    GAW.grand_avg_inhibSDc          =  GAW.grand_avg_inhibSDc       / total_alpha_cycles; 
    GAW.grand_avg_FRSDc             =  GAW.grand_avg_FRSDc          / total_alpha_cycles; 
    GAW.grand_avg_ExternSDc         =  GAW.grand_avg_ExternSDc      / total_alpha_cycles; 
    GAW.grand_avg_ExternInhSDc		=  GAW.grand_avg_ExternInhSDc	/ total_alpha_cycles; 
    GAW.grand_avg_RecurrentSDc 		=  GAW.grand_avg_RecurrentSDc 	/ total_alpha_cycles; 
    GAW.grand_avg_inputSDc  		=  GAW.grand_avg_inputSDc       / total_alpha_cycles; 
    GAW.grand_avg_inputInhSDc		=  GAW.grand_avg_inputInhSDc	/ total_alpha_cycles; 
    GAW.grand_avg_ExcSynSDc 		=  GAW.grand_avg_ExcSynSDc      / total_alpha_cycles; 
    GAW.grand_avg_GlobalSDc 		=  GAW.grand_avg_GlobalSDc      / total_alpha_cycles; 
    GAW.grand_avg_InhSynSDc 		=  GAW.grand_avg_InhSynSDc      / total_alpha_cycles; 
    GAW.grand_avg_FR_inhSDc 		=  GAW.grand_avg_FR_inhSDc      / total_alpha_cycles; 
    GAW.grand_avg_BOLDexcSDc		=  GAW.grand_avg_BOLDexcSDc     / total_alpha_cycles; 
    GAW.bin1                        =  GAW.bin1                     / lenb1; 
    GAW.bin2                        =  GAW.bin2                     / lenb2; 
    GAW.bin3                        =  GAW.bin3                     / lenb3; 
    GAW.bin4                        =  GAW.bin4                     / lenb4; 
    GAW.bin5                        =  GAW.bin5                     / lenb5; 
    GAW.bin6                        =  GAW.bin6                     / lenb6; 


    % Iterate again over all GAWnode files to compute squared deviation
    % from mean
    for ii = 1:length(GAWnode_filenames)

        % Load GAWnode file
        load(GAWnode_filenames{ii});

        % Compute squared deviation from mean
        SD_GAW.grand_avg_alphaSDc		= SD_GAW.grand_avg_alphaSDc     + sum(( GAWnode.grand_avg_alphaSDc      - GAW.grand_avg_alphaSDc 	).^2);
        SD_GAW.grand_avg_inhibSDc 		= SD_GAW.grand_avg_inhibSDc 	+ sum(( GAWnode.grand_avg_inhibSDc      - GAW.grand_avg_inhibSDc 	).^2);
        SD_GAW.grand_avg_FRSDc			= SD_GAW.grand_avg_FRSDc        + sum(( GAWnode.grand_avg_FRSDc         - GAW.grand_avg_FRSDc	 	).^2);
        SD_GAW.grand_avg_ExternSDc		= SD_GAW.grand_avg_ExternSDc	+ sum(( GAWnode.grand_avg_ExternSDc     - GAW.grand_avg_ExternSDc	).^2);
        SD_GAW.grand_avg_ExternInhSDc	= SD_GAW.grand_avg_ExternInhSDc + sum(( GAWnode.grand_avg_ExternInhSDc  - GAW.grand_avg_ExternInhSDc).^2);
        SD_GAW.grand_avg_RecurrentSDc 	= SD_GAW.grand_avg_RecurrentSDc + sum(( GAWnode.grand_avg_RecurrentSDc  - GAW.grand_avg_RecurrentSDc).^2);
        SD_GAW.grand_avg_inputSDc  		= SD_GAW.grand_avg_inputSDc  	+ sum(( GAWnode.grand_avg_inputSDc      - GAW.grand_avg_inputSDc  	).^2);
        SD_GAW.grand_avg_inputInhSDc	= SD_GAW.grand_avg_inputInhSDc	+ sum(( GAWnode.grand_avg_inputInhSDc   - GAW.grand_avg_inputInhSDc	).^2);
        SD_GAW.grand_avg_ExcSynSDc 		= SD_GAW.grand_avg_ExcSynSDc 	+ sum(( GAWnode.grand_avg_ExcSynSDc     - GAW.grand_avg_ExcSynSDc 	).^2);
        SD_GAW.grand_avg_GlobalSDc 		= SD_GAW.grand_avg_GlobalSDc 	+ sum(( GAWnode.grand_avg_GlobalSDc     - GAW.grand_avg_GlobalSDc 	).^2);
        SD_GAW.grand_avg_InhSynSDc 		= SD_GAW.grand_avg_InhSynSDc 	+ sum(( GAWnode.grand_avg_InhSynSDc     - GAW.grand_avg_InhSynSDc 	).^2);
        SD_GAW.grand_avg_FR_inhSDc 		= SD_GAW.grand_avg_FR_inhSDc 	+ sum(( GAWnode.grand_avg_FR_inhSDc     - GAW.grand_avg_FR_inhSDc 	).^2);
        SD_GAW.grand_avg_BOLDexcSDc		= SD_GAW.grand_avg_BOLDexcSDc	+ sum(( GAWnode.grand_avg_BOLDexcSDc    - GAW.grand_avg_BOLDexcSDc	).^2);
        SD_GAW.bin1                     = SD_GAW.bin1                   + sum(( GAWnode.bin1                    - GAW.bin1                  ).^2);
        SD_GAW.bin2                     = SD_GAW.bin2                   + sum(( GAWnode.bin2                    - GAW.bin2                  ).^2);
        SD_GAW.bin3                     = SD_GAW.bin3                   + sum(( GAWnode.bin3                    - GAW.bin3                  ).^2);
        SD_GAW.bin4                     = SD_GAW.bin4                   + sum(( GAWnode.bin4                    - GAW.bin4                  ).^2);
        SD_GAW.bin5                     = SD_GAW.bin5                   + sum(( GAWnode.bin5                    - GAW.bin5                  ).^2);
        SD_GAW.bin6                     = SD_GAW.bin6                   + sum(( GAWnode.bin6                    - GAW.bin6                  ).^2);

    end

    % Compute standard deviation
    SD_GAW.grand_avg_alphaSDc           = sqrt(SD_GAW.grand_avg_alphaSDc	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_inhibSDc           = sqrt(SD_GAW.grand_avg_inhibSDc 	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_FRSDc              = sqrt(SD_GAW.grand_avg_FRSDc       / (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_ExternSDc          = sqrt(SD_GAW.grand_avg_ExternSDc	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_ExternInhSDc		= sqrt(SD_GAW.grand_avg_ExternInhSDc/ (total_alpha_cycles - 1) );
    SD_GAW.grand_avg_RecurrentSDc       = sqrt(SD_GAW.grand_avg_RecurrentSDc/ (total_alpha_cycles - 1) );
    SD_GAW.grand_avg_inputSDc           = sqrt(SD_GAW.grand_avg_inputSDc  	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_inputInhSDc		= sqrt(SD_GAW.grand_avg_inputInhSDc	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_ExcSynSDc          = sqrt(SD_GAW.grand_avg_ExcSynSDc 	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_GlobalSDc          = sqrt(SD_GAW.grand_avg_GlobalSDc 	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_InhSynSDc          = sqrt(SD_GAW.grand_avg_InhSynSDc 	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_FR_inhSDc          = sqrt(SD_GAW.grand_avg_FR_inhSDc 	/ (total_alpha_cycles - 1) ); 
    SD_GAW.grand_avg_BOLDexcSDc         = sqrt(SD_GAW.grand_avg_BOLDexcSDc	/ (total_alpha_cycles - 1) ); 
    SD_GAW.bin1                         = sqrt(SD_GAW.bin1                  / (total_alpha_cycles - 1) ); 
    SD_GAW.bin2                         = sqrt(SD_GAW.bin2                  / (total_alpha_cycles - 1) ); 
    SD_GAW.bin3                         = sqrt(SD_GAW.bin3                  / (total_alpha_cycles - 1) ); 
    SD_GAW.bin4                         = sqrt(SD_GAW.bin4                  / (total_alpha_cycles - 1) ); 
    SD_GAW.bin5                         = sqrt(SD_GAW.bin5                  / (total_alpha_cycles - 1) ); 
    SD_GAW.bin6                         = sqrt(SD_GAW.bin6                  / (total_alpha_cycles - 1) ); 


    % Compute standard error of the mean
    GAW.inhibSD               =    SD_GAW.grand_avg_inhibSDc        /     sqrt(total_alpha_cycles);
    GAW.RecurrentSD           =    SD_GAW.grand_avg_RecurrentSDc    /     sqrt(total_alpha_cycles);
    GAW.GlobalSD              =    SD_GAW.grand_avg_GlobalSDc       /     sqrt(total_alpha_cycles);
    GAW.fmriexcSD             =    SD_GAW.grand_avg_BOLDexcSDc      /     sqrt(total_alpha_cycles);
    GAW.alphaSD               =    SD_GAW.grand_avg_alphaSDc        /     sqrt(total_alpha_cycles);
    GAW.excSynSD              =    SD_GAW.grand_avg_ExcSynSDc       /     sqrt(total_alpha_cycles);
    GAW.inhSynSD              =    SD_GAW.grand_avg_InhSynSDc       /     sqrt(total_alpha_cycles);
    GAW.extSD                 =    SD_GAW.grand_avg_ExternSDc       /     sqrt(total_alpha_cycles);
    GAW.extinhSD              =    SD_GAW.grand_avg_ExternInhSDc    /     sqrt(total_alpha_cycles);
    GAW.inhSD                 =    SD_GAW.grand_avg_InhSynSDc       /     sqrt(total_alpha_cycles);
    GAW.excSD                 =    SD_GAW.grand_avg_RecurrentSDc    /     sqrt(total_alpha_cycles);
    GAW.inpExcSD              =    SD_GAW.grand_avg_inputSDc        /     sqrt(total_alpha_cycles);
    GAW.inpInhSD              =    SD_GAW.grand_avg_inputInhSDc     /     sqrt(total_alpha_cycles)/1.4;
    GAW.FRExcSD               =    SD_GAW.grand_avg_FRSDc           /     sqrt(total_alpha_cycles);
    GAW.FRInhSD               =    SD_GAW.grand_avg_FR_inhSDc       /     sqrt(total_alpha_cycles);
    GAW.bin1SD                =    SD_GAW.bin1                      /     sqrt(lenb1);
    GAW.bin2SD                =    SD_GAW.bin2                      /     sqrt(lenb2);
    GAW.bin3SD                =    SD_GAW.bin3                      /     sqrt(lenb3);
    GAW.bin4SD                =    SD_GAW.bin4                      /     sqrt(lenb4);
    GAW.bin5SD                =    SD_GAW.bin5                      /     sqrt(lenb5);
    GAW.bin6SD                =    SD_GAW.bin6                      /     sqrt(lenb6);

end
