function average_perc_difference_replicates=compute_percent_difference_between_replicates(data_output)

%Computes the average percent difference between all replicates in the data
%using the diff_measurements./mean_measurements
counter=1;

for iData_Output=1:length(data_output)
    
    %Make sure there is more than one replicate
    if length(data_output(iData_Output).values)>1
        
        Query_Strain_vals=data_output(iData_Output).values;
        
        %% x and y values for each of the strains
        x=Query_Strain_vals(1);
        y=Query_Strain_vals(2);
        
        Query_Strain_Difference=abs(x-y);
        
        Query_Strain_Average=abs((x+y)./2);
        %Query_Strain_Average=(x+y)./2;
        
        Query_percent(counter)=Query_Strain_Difference./Query_Strain_Average;
        
        counter=counter+1;
        
    end
    
end

average_perc_difference_replicates=nanmean(Query_percent);
