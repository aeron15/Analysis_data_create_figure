
%% Compute the 10% above the mean and 10% below the mean

load('data/setpoints_normalized');

%Compute the 95% confidence interval for estimate

BC187_vals_vector=cell2mat(setpoints_normalized(:,1));

BC187_mean=nanmean(BC187_vals_vector);

BC187_std=nanstd(BC187_vals_vector);


% Standard deviation and length of the vector
%SEM = std(x)/sqrt(length(x));  
BC187_len=sum(~isnan(BC187_vals_vector));

BC187_SEM = BC187_std/sqrt(BC187_len); 


higher_bound= 1.1* BC187_mean;
lower_bound=0.9 * BC187_mean;
vline(higher_bound)
vline(lower_bound)
title ('BC187 measurements across replicates')

idx_to_remove=find(~(higher_bound * BC187_mean & BC187_vals_vector < lower_bound));
setpoints_normalized(idx_to_remove,:)=[];