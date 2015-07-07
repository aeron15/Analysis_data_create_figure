function driver_main_figures

%DRIVER_MAIN_FIGURES plots the main figures of the paper in a nice
%format for publication

load ../outputFigures/data_output_figure_1.mat
plot_figure_mean_error_bar_side_histogram(data_output,'file_append','figure1')

%Correlation genetic distance and set point of induction using RAD-seq data
NaturalIsolates_correlation=compute_correlation_genetic_distance_set_point_induction(data_output,loc);

%%
load ../outputFigures/data_output_figure_3.mat
plot_figure_mean_error_bar(data_output,'file_append','figure3','figure_size','AlleleReplacement')

%%
load ../outputFigures/data_output_figure_4.mat
plot_figure_mean_error_bar(data_output,'file_append','figure4')

%Compute correlation coefficient between natural isolates and allele
%replacements
[NaturalIsolatesSwaps_CorrelationCoefficient,NaturalIsolatesSwaps_PValueCorrelation]=compute_correlation_natural_isolates_allele_replacements;

%% Make GAL3 hemizygous hybrid

load ../outputFigures/data_output_figure_GAL3HH.mat
plot_figure_mean_error_bar(data_output,'file_append','figureGAL3HH','figure_size','HH')

%% Make SOK1 hemizygous hybrid

load ../outputFigures/data_output_figure_SOK1HH.mat
plot_figure_mean_error_bar(data_output,'file_append','figureSOK1HH','figure_size','HH')


%% Make figure 4
close all;
load('all_strains_names')
load('all_strains_vals_vector')

fig4(all_strains_vals_vector, all_strains_names);

%% Make figure 5  

close all;
load('all_strains_names')
load('all_strains_vals_vector')
fig5(all_strains_vals_vector, all_strains_names);

