function [vector,strains]=compute_mean_data_output(data_output)

%COMPUTE_MEAN_DATA_OUTPUT of each all strains and exports the strain names

strains={data_output.strain};

for iStrain = 1:length(strains)
    
    QueryStrain_vals=data_output(iStrain).values;
    QueryStrain_mean=mean(QueryStrain_vals);
    QueryStrain_std=std(QueryStrain_vals);
    
    
    data_output_mean(iStrain).strain=data_output(iStrain).strain;
    data_output_mean(iStrain).values=QueryStrain_mean;
    
    vector(iStrain)=QueryStrain_mean;    
    
end


end