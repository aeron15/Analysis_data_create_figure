function [converted_names] = natural_isos_ref(input_names)

% The strain_name_conv table takes all of the strains that are normally
% processed in the lab from a 3 column excel file and converts them to a
% usable form.

% Most strains have names that were appended to them during strain
% construction (column 1), most of these names have the same first
% characters so they can also be written as a wildcard (column 2), and in
% order to be processed easily (column 3) contains name descriptions that
% can be used for labeling.

% To update strain list, find the original 'strain_conv_table.xlsx' file.

% Created by KL 20150227

%% LOAD strain conversion table

load('Strain_conversion_table.mat');
master_list = output;

%% Determine if input_names isacell

TF = iscell(input_names);

%% Assign label name from construct or wildcard name

construct_name = master_list(:,1);
wildcard_name = master_list(:,2);
label_name = master_list(:,3);
natural_iso_ref = master_list(:,4);

for iNum = 1:length(input_names)
    
    if TF == 1

        new_name = strcmp(input_names{iNum}, construct_name);
        idx = find(new_name(:,1)==1);

        try
            temp_label{iNum} = natural_iso_ref{idx};
        end

        try
            new_name = strcmp(input_names{iNum}, wildcard_name);
            idx_2 = find(new_name(:,1)==1);

            temp_label{iNum} = natural_iso_ref{idx_2};
        catch
            temp_label{iNum} = [];
        end
        
        converted_names = temp_label;
        
    else
   
         new_name = strcmp(input_names, construct_name);
         idx = find(new_name(:,1)==1);

            try
                temp_label = label_name{idx};
            catch
                new_name = strcmp(input_names, wildcard_name);
                idx_2 = find(new_name(:,1)==1);

                temp_label = label_name{idx_2};
            end

            converted_names = temp_label;
            
    end
    
end
   
    


