function average_standard_deviation=compute_average_standard_deviation(data_output)

%Computes the average standard deviation in the data
%using the diff_measurements./mean_measurements

counter=1;

for iData_Output=1:length(data_output)
    
    Query_Strain_vals=data_output(iData_Output).values;
    standard_deviation_vector(counter)=nanstd(data_output(iData_Output).values);
    counter=counter+1;
    
end

average_standard_deviation=nanstd(standard_deviation_vector);

end


