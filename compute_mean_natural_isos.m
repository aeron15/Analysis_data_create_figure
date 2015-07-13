function [GLUmean_naturalIsos,medianOfDist,GAL_GLU_ratio]=compute_mean_natural_isos(data_output)

%Computes the mean of the distribution of natural isolates and other
%descriptive features of the distribution

counter=1;

for iData_Output=1:length(data_output)
    
    QueryStrain_vals=data_output(iData_Output).values;
    
    QuerySetpoint_mean(counter) = nanmean(QueryStrain_vals);
    QuerySetpoint_std(counter) = nanstd(QueryStrain_vals);
    
    counter=counter+1;
     
end

%% Describing the distrubution

meanOfDist = nanmean(QuerySetpoint_mean);
medianOfDist = median(QuerySetpoint_mean);
stdOfDist = nanstd(QuerySetpoint_mean);

%Convert mean to glucose space
GLUmean_naturalIsos = 2^meanOfDist;
GAL_GLU_ratio = .125/GLUmean_naturalIsos;

[distribution_min, loc1] = min(QuerySetpoint_mean);
dist_min_GLU = 2^distribution_min;
error_std_min = QuerySetpoint_std(loc1);

[distribution_max, loc2] = max(QuerySetpoint_mean);
dist_max_GLU = 2^distribution_max;
error_std_max = QuerySetpoint_std(loc2);

end

