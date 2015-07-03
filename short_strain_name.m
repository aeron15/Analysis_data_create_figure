function converted_names = short_strain_name(input_names)

% The SHORT_STRAIN_NAME converts the strains that have different '_'
% appendages to the same unique strain name identifier.

% To update strain list, find the original 'short_strain_conv.xlsx' file.

% Created by KL 20150703


%% LOAD strain conversion table

load('unique_strains_figure_table.mat');
master_list = output;

%% Determine if input_names isacell

TF = iscell(input_names);

%% Assign label name from construct or wildcard name

construct_name = master_list(:,1);
other_name = master_list(:,2);

if TF==1

    for iNum = 1:length(input_names)

        new_name = strcmp(input_names{iNum}, construct_name);
        idx = find(new_name(:,1)==1);

            try
                temp_label{iNum} = other_name{idx};
            catch
                temp_label{iNum} = input_names{iNum};
            end

    end
    
    converted_names = temp_label;
else
    
     new_name = strcmp(input_names, construct_name);
     idx = find(new_name(:,1)==1);

        try
            temp_label = other_name{idx};
        catch
            temp_label = input_names;
        end
        
        converted_names = temp_label;
        
end
    
    

    
    


