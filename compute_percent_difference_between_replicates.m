function average_perc_difference_replicates=compute_percent_difference_between_replicates(data_output)

%COMPUTE_PERCENT_DIFFERENCE_BETWEEN_REPLICATES
%Computes the average percent difference between all replicates in the data
%using the formula diff_measurements./mean_measurements

counter=1;

for iData_Output=1:length(data_output)
    
    %Make sure there is more than one replicate
    if length(data_output(iData_Output).values)>1
        
        QueryStrain_vals=data_output(iData_Output).values;
        
        %% x and y values for each of the strains
        x=QueryStrain_vals(1);
        y=QueryStrain_vals(2);
        
        QueryStrain_Difference=abs(x-y);
        
        QueryStrain_Average=abs((x+y)./2);
        %Query_Strain_Average=(x+y)./2;
        
        QueryStrain_percent(counter)=QueryStrain_Difference./QueryStrain_Average;
        
        counter=counter+1;
        
    end
    
end

average_perc_difference_replicates=nanmean(QueryStrain_percent);
