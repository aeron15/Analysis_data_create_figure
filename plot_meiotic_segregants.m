function plot_meiotic_segregants()

load data_meiotic_segregants.mat

%%
plot_distribution(meioticSegregants_mean)
filename=['Distribution_meiotic_segregants_mean.pdf'];
export_fig_specific_path(filename,'-pdf',  '-transparent', '-nocrop')

%%
plot_distribution(meioticSegregants_area)
filename=['Distribution_meiotic_segregants_area.pdf'];
export_fig_specific_path(filename,'-pdf',  '-transparent', '-nocrop')


%% The area metric is much more biased
figure;
scatterhist(meioticSegregants_mean,meioticSegregants_area)
axis square;
filename=['Distribution of meiotic segregants_mean_vs_area_more_events.pdf'];
export_fig_specific_path(filename,'-pdf',  '-transparent', '-nocrop')
