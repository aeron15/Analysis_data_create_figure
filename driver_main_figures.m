function driver_main_figures

%DRIVER_MAIN_FIGURES plots the main figures of the paper in a nice
%format for publication

load data_output_figure_1.mat
plot_figure_mean_error_bar_side_histogram(data_output,'file_append','figure1')

%%
load data_output_figure_3.mat
plot_figure_mean_error_bar(data_output,'file_append','figure3','figure_size','AlleleReplacement')

%%
load data_output_figure_4.mat
plot_figure_mean_error_bar(data_output,'file_append','figure4')

%% Print figure 6
load('all_strains_names')
load('all_strains_vals_vector')

fig5(all_strains_vals_vector, all_strains_names);
fig6(all_strains_vals_vector, all_strains_names);
%% Print GAL3 hemizygous hybrid

load data_output_figure_GAL3HH.mat
plot_figure_mean_error_bar(data_output,'file_append','figureGAL3HH','figure_size','HH')

%% Print SOK1 hemizygous hybrid

load data_output_figure_SOK1HH.mat
plot_figure_mean_error_bar(data_output,'file_append','figureSOK1HH','figure_size','HH')

%%

