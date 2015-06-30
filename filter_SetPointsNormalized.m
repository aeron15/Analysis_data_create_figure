function setpoints_normalized = filter_SetPointsNormalized(setpoints_normalized)

%FILTER_SETPOINTSNORMALIZED filters data using the values of BC187 and
%plots the filtered data

% Compute standard deviation and length of the vector

BC187_vals_vector=cell2mat(setpoints_normalized(:,1));

BC187_mean=nanmean(BC187_vals_vector);

BC187_std=nanstd(BC187_vals_vector);


%% Compute the 95% confidence interval for estimate (this part is actually not used for any filtering

BC187_len=sum(~isnan(BC187_vals_vector));

BC187_SEM = BC187_std/sqrt(BC187_len);  %SEM = std(x)/sqrt(length(x));  

ts = tinv([0.01  0.99],BC187_len-1);      % T-Score

CI = BC187_mean + ts*BC187_SEM;

%% Compute the 10% above the mean and 10% below the mean and remove bad data

higher_bound= 1.1* BC187_mean;
lower_bound=0.89 * BC187_mean;

idx_to_remove=~(higher_bound * BC187_mean & BC187_vals_vector < lower_bound);
setpoints_normalized(idx_to_remove,:)=[];
%% Plot filtered data (Check if the plot is done for galactose or glucose titrations)

%plot_hist_BC187_vals(BC187_vals_vector,higher_bound,lower_bound)

end