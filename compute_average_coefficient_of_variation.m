function AverageCoefficientVariation=compute_average_coefficient_of_variation(data_output)

%Computes the average standard deviation in the data
%using the diff_measurements./mean_measurements

counter=1;

for iData_Output=1:length(data_output)
    
    QueryStrain_vals=data_output(iData_Output).values;
    %Computed weighted coefficient of variation
    
    CoefficientVariation_vector(counter)=nanstd(QueryStrain_vals)./abs(nanmean(QueryStrain_vals));
    Weights(counter)=length(QueryStrain_vals);
    
    counter=counter+1;
    
    
end

%Remove cases where the variation is zero those are cases where there is
%only one replicate
AllStrains=({data_output.strain});
StrainsOneReplicate_idx=find(CoefficientVariation_vector==0);
AllStrains(StrainsOneReplicate_idx);

%Remove strains with 1 replicate and NaNs
strainsOneReplicate_idx=(CoefficientVariation_vector==0 | isnan(CoefficientVariation_vector));

CoefficientVariation_vector(strainsOneReplicate_idx)=[];
Weights(strainsOneReplicate_idx)=[];

%Weighted coefficient of variation
AverageCoefficientVariation=wmean(CoefficientVariation_vector,Weights);

end
