function plot_hist_BC187_vals(BC187_vals_vector,higher_bound,lower_bound)

%PLOT_HIST_BC187_VALS plots all the values of BC187 to show the data that
%is saved and the data that is discarded

BC187_mean=nanmean(BC187_vals_vector);
 
hfig=figure;
hist(BC187_vals_vector,20)
vline(BC187_mean,'green')
vline(higher_bound)
vline(lower_bound)
title ('BC187 measurements across replicates')

Set_fig_RE(hfig,9,9,20)
xlabel('Set point of induction of BC187')


filename='BC187_vals_histogram';
export_fig(filename, '-pdf','-transparent','-nocrop');
