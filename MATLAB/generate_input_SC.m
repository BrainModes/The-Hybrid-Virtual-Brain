% Converts MATLAB format SC files into a txt file format that can be read 
% by the hybrid brain network model implementation
%
% USAGE: 
% generate_input_SC(SC, output_filestem)
%
% INPUTS:
% SC                - [N,N] array that contains connection weights for all 
%                     NxN region pairs
% output_filestem   - String that becomes the prefix of the three output
%                     files
%                 
% OUTPUTS:
% The hybrid brain network model implementation needs three input text
% files. All three filenames must contain identical prefixes and this 
% prefix must be supplied as an argument for the brain model.
% Suffix _SC_strengths.txt: contains connection weights,
% suffix _SC_distances.txt: contains connection distances,
% suffix _SC_regionids.txt: contains the region numbers of the source and
% target regions for each region of the brain model


function generate_input_SC(SC, output_filestem)

    % Normalize SC by dividing through maximum value
    SC          = SC ./ max(SC(:));
    SCsize      = size(SC,1);
    
    % OPTIONAL: Set low coupling values to zero
    %sumSC       = sum(SC);
    %for ii = 1:SCsize
    %    SC(ii, find(SC(ii,:) < sumSC(ii)*0.05) ) = 0;
    %end
        
    % Generate filenames
    sc_cap_file  = [output_filestem '_SC_strengths.txt'];
    sc_dist_file = [output_filestem '_SC_distances.txt'];
    sc_id_file   = [output_filestem '_SC_regionids.txt'];
    
    % Write number of nodes as header line into SC files
    dlmwrite(sc_cap_file,   SCsize);
    dlmwrite(sc_dist_file,  SCsize);
    dlmwrite(sc_id_file,    SCsize);

    % Write maximum distance as 2nd header line into dist file
    % Time-delays are turned off, max-distance is set to 1
    maxdist     =  1;
    dlmwrite(sc_dist_file,maxdist,'delimiter',' ','-append');

    % Write connection weights, distances and input-region IDs
    for ii = 1:SCsize,
        % Select all non-zero connections
        inpregs     =   find(SC(ii,:)>0);
        inpcaps     =   SC(ii,inpregs);
        inpdists    =   sc_dist(ii,inpregs);

        % Write connectivity description for each node
        % Convert from Matlab-style 1-based numbering to C-style 0-based 
        % numbering
        inpregs     =   inpregs-1; 
        cap_line    =   [(ii-1) length(inpregs)];
        dist_line   =   [(ii-1) length(inpregs)];
        inp_line    =   [(ii-1) length(inpregs)];
        dlmwrite(sc_cap_file,   cap_line,'delimiter',' ','-append');
        dlmwrite(sc_dist_file,  dist_line,'delimiter',' ','-append');
        dlmwrite(sc_id_file,    inp_line,'delimiter',' ','-append');

        % Write actual connectivity information
        dlmwrite(sc_cap_file,   inpcaps,'delimiter',' ','-append','precision','%.8f');
        dlmwrite(sc_dist_file,  inpdists,'delimiter',' ','-append','precision','%.8f');
        dlmwrite(sc_id_file,    inpregs,'delimiter',' ','-append');
    end

end
