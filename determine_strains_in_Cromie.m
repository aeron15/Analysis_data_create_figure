function QueryStrains_counter=determine_strains_in_Cromie()

%DETERMINE_STRAINS_IN_CROMIE determines the number of strains that we use
%in the paper that were also sequenced in the Cromie paper. 
%Determine which strains are missing

load('data_output_figure_1')
load('Cromie_distances_2.mat')


Strain_names={data_output.strain};
QueryStrains_counter=0;

for iStrain=1:length(Strain_names)
    
    QueryStrain=Strain_names{iStrain};
    
    %Check if replacing is needed

    QueryStrain=replace_strain(QueryStrain);

    %Find index of the strain
    
    CurrStrain = regexp(CromieStrains_names, QueryStrain);
    cs = cellfun(@isempty,CurrStrain);
    QueryStrain_idx = find(cs==0);
    
    if QueryStrain_idx
        
        QueryStrains_counter=QueryStrains_counter+1;
        
    else
        
        QueryStrain
        
    end
    
end




