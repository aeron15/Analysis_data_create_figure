function compute_setpoints_reference_BC187()

%COMPUTE_SETPOINTS_REFERENCE_BC187

%Input should be the set points normalized in data
% To plot data look a the script driver_main_figures.m

%% Load data

load(['../data/setpoints_normalized_area.mat']);

%% Filter data using 2 standard deviations from the median reference BC187 value
setpoints_normalized=filterData_onlyforNaturalIsolates(setpoints_normalized);

%% Create variable equivalent to dif_sp.mat from the plot_all_figs_1

all_strains_vals_vector=cell2mat(setpoints_normalized(:,2));
all_strains_names=setpoints_normalized(:,3);

save('../outputFigures/all_strains_vals_vector','all_strains_vals_vector');
save('../outputFigures/all_strains_names','all_strains_names');

%% Setpoints of all strains in study

all_strains  = {'Y55*'; 'NCYC110*'; 'L_1528*'; 'DBVPG6044*';
    'Y12_SGRP*'; 'W303*'; 'i378604X*'; 'DBVPG1373*';
    'YIIc17_E5*'; 'UWOPS87_2421*'; 'YPS163*'; 'CLIB215*';
    'CLIB324*'; 'NC_02*'; 'PW5*'; 'YS4*';
    'T7*'; 'Y9_WashU*'; 'UWOPS03_4614*'; 'IL_01*';
    'M22*'; 'DBVPG6765*'; 'YPS128*'; 'DBVPG1788*';
    'DBVPG1853*'; 'L_1374*'; 'DBVPG1106*'; 'YJM421*';
    'Bb32*'; 'YJM428*'; 'UWOPS05_2272*'; 'DBVPG6040*';
    'YJM653*'; 'UC5*'; 'YPS1009*'; 'CLIB382*';
    'WE372*'; 'YJM975*'; 'I_14*'; 'YJM981*';
    'Y12_WashU*'; 'FL100*'; 'i273614N*'; 'YPS606*';
    'BC187*'; 'YJM978*'; 'S288C*';'RY16*'; 'RYB53*'; 'RYB59*'; 'RYB65*'; 'RYB66*'; 'RYB28*';
    'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';
    'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*';
    'RYC45*';'RYC58*';'RYC49*'; 'RYC50*';'RYC51*'; 'RYC59_1*';'RYC52*';'RYC60*';'RYC62*'; 'RYB92*'; 'RYC72*';
    'RYD25*'; 'RYD27*'; 'RYD28*'; 'RYD30*'; 'RYD31*'; 'RYB59*'; 'RYB53*'};

rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);
filename='All_data';

[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

%Compute coefficient of variation across all the paper
AllData_CoefficientOfVariation=compute_average_coefficient_of_variation(data_output);

save('../outputFigures/data_output_figure_glucose_titration','data_output');




%% Figure 1. Setpoints of natural isolates

%Removed CLIB215
all_strains  = {'Y55*'; 'NCYC110*'; 'L_1528*'; 'DBVPG6044*';
    'Y12_SGRP*'; 'W303*'; 'i378604X*'; 'DBVPG1373*';
    'YIIc17_E5*'; 'UWOPS87_2421*'; 'YPS163*';
    'CLIB324*'; 'NC_02*'; 'PW5*'; 'YS4*';
    'T7*'; 'Y9_WashU*'; 'UWOPS03_4614*'; 'IL_01*';
    'M22*'; 'DBVPG6765*'; 'YPS128*';
    'DBVPG1853*'; 'L_1374*'; 'DBVPG1106*'; 'YJM421*';
    'Bb32*'; 'YJM428*'; 'UWOPS05_2272*'; 'DBVPG6040*';
    'YJM653*'; 'UC5*'; 'YPS1009*'; 'CLIB382*';
    'WE372*'; 'YJM975*'; 'I_14*'; 'YJM981*';
    'Y12_WashU*'; 'FL100*'; 'i273614N*'; 'YPS606*';
    'BC187*'; 'YJM978*'; 'S288C*'};

rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);
filename='Fig_1_natural_isolates';

[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('../outputFigures/data_output_figure_1','data_output','loc');
save('../outputFigures/data_output_natural_isolates_glucose_titration','data_output','loc');

%Compute coefficient of variation across all the paper
NaturalIsolates_CoefficientOfVariation=compute_average_coefficient_of_variation(data_output);


%>>>>>LOG ENTRIES
add_entry_log('Number of natural isolates',length({data_output.strain}));
add_entry_log('Natural Isolates coefficient of variation across replicates', NaturalIsolates_CoefficientOfVariation);

%Correlation genetic distance and set point of induction using RAD-seq data
NaturalIsolates_correlation=compute_correlation_genetic_distance_set_point_induction(data_output,loc);
add_entry_log('Correlation between genetic distance and set point difference', NaturalIsolates_correlation);

%Highlight specific points of the data in the correlation plot as a
%reference
%plot_correlation_highlight_points()

%QueryStrains_counter=determine_strains_in_Cromie();

% Compute difference between BC187 and YJM978
strain1='BC187'; strain2='YJM978';
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2);

%Range of the set points on figure 1 (Natural Isolates)
strain1= data_output(loc(1)).strain;
strain2= data_output(loc(end)).strain;
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2);

% Average standard deviation of natural isolates set point of induction
NaturalIsolates_AverageStandardDeviation=compute_average_standard_deviation(data_output);

% Average standard deviation of natural isolates set point of induction
NaturalIsolates_CoefficientOfVariation=compute_average_coefficient_of_variation(data_output);

%Range of the natural isolates strains on figure 4
strain1='YJM421';
strain2='DBVPG1373';
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2);

%Range of the strains on figure 5
strain1='YJM978';
strain2='I-14';
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2);


%% Figure 3. Allele swaps BC187-YJM978
%strains = {'RY16*', 'RYB53*', 'RYB59*', 'RYB65*', 'RYB66*', 'RYB28*'};
strains = {'RYB53*', 'RYB59*', 'RYB65*', 'RYB66*'};

filename='Fig_3_allele_swaps';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('data_output_figure_3','data_output');

%compute_conversion_rate

InterStrain_distance=abs(mean(data_output(2).values)-mean(data_output(3).values));

%data_output(1).values=data_output(1).values([1:3 5]);

YJM978background_differences=abs(mean(data_output(2).values)-mean(data_output(1).values))
BC187background_differences=abs(mean(data_output(3).values)-mean(data_output(4).values));

%Compute percent conversion strains
YJM978bg_PercentConversion=YJM978background_differences./InterStrain_distance;
BC187bg_PercentConversion=BC187background_differences./InterStrain_distance;
add_entry_log('Conversion in the YJM978 background', NaturalIsolates_correlation);
add_entry_log('Conversion in the YJM978 background', NaturalIsolates_correlation);

%Range of the strains on figure 3
strain1='GAL3-BC (YJM978)';
strain2='GAL3-YJM (YJM978)';
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2);

%Range of the strains on figure 3
strain1='GAL3-BC (BC187) ';
strain2='GAL3-YJM (BC187) ';
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2);

%% Figure 3 test heterologos locus effect
strains = {'RY16*', 'RYB53*', 'RYB59*', 'RYB65*', 'RYB66*', 'RYB28*'};

filename='Fig_3_allele_swaps';
[fig3,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

BC_het_effect=mean([fig3(4).values;fig3(6).values]);
YJ_het_effect=mean([fig3(1).values;fig3(3).values]);

Delta_GAL3=BC_het_effect-YJ_het_effect;

set_1=[fig3(1).values;fig3(6).values-Delta_GAL3];%haploids
set_2=[fig3(3).values;fig3(4).values-Delta_GAL3];

%T-TEST
[h,p,ci,stats]=ttest2(set_1,set_2);

%>>>>LOG ENTRY
add_entry_log('Heterologous locus effect T-test result', p);

%% Figure 4. Natural Isolate ORF swaps into YJM978

%strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*','RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*'};
strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*', 'RYD25*', 'RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*'};
%Remove RYB92 which is the S288C allele
filename='Fig_4_YJ_bg_Diff_alleles';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('../outputFigures/data_output_figure_4','data_output');

%Range of variation of the allele replacements
strain1= data_output(loc(1)).strain;
strain2= data_output(loc(end)).strain;
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2)

%>>>>LOG ENTRY
add_entry_log('Range variation Allele Swaps Figure 4 background YJM978', FoldDifferenceMean);
add_entry_log('Error range of variation fold Difference', ErrorFoldDifference);



%Compute correlation coefficient between natural isolates and allele
%replacements
[NaturalIsolatesSwaps_CorrelationCoefficient,NaturalIsolatesSwaps_PValueCorrelation]=compute_correlation_natural_isolates_allele_replacements;
add_entry_log('Correlation natural isoaltes and allele swaps', ErrorFoldDifference);
add_entry_log('P-value Correlation natural isoaltes and allele swaps', NaturalIsolatesSwaps_PValueCorrelation);


%% Figure 5 BC187 alleles
%strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'; 'RYD50*'; 'RYD53*'};
strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'};

filename='Fig5_BC';
[data_output,loc]=make_dot_plot(strains_BC, all_strains_vals_vector, all_strains_names, filename);
StrainsWithBC187Allele_names={data_output.strain}';

%Range of variation of the allele replacements
% strain1= 'GAL3-BC (YJM978)';
% strain2= 'GAL3-BC187 (I-14)';
strain1= data_output(loc(1)).strain;
strain2= data_output(loc(end)).strain;
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2)


%% Figure 5 YJM978 alleles
%strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*'};
strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'};
filename='Fig5_YJ_alelle';
[data_output,loc]=make_dot_plot(strains_YJM, all_strains_vals_vector, all_strains_names, filename);
StrainsWithYJM978Allele_names={data_output.strain}';

strain1= data_output(loc(1)).strain;
strain2= data_output(loc(end)).strain;
[~,~,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2)


%% Figure 5 all allele replacements of YJM978 and BC187
close all;
%strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*'};
strains = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'};

filename='Fig5_BC_YJ';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('../outputFigures/data_output_figure_5','data_output');

%Determine fold differences between pairs of strains
[AlleleReplacementBackgrounds_mean,AlleleReplacementBackgrounds_std,AlleleReplacementBackgrounds_sem]=compute_fold_difference_across_backgrounds(data_output,StrainsWithBC187Allele_names,StrainsWithYJM978Allele_names)

%% SI figures. Hemizygous hybrid strains
close all;
strains={'RYB22';
    'RYB23';
    'RYB24'}
%'RYB42'}
filename='GAL3_HH';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('../outputFigures/data_output_figure_GAL3HH','data_output');

[GAL3YJHH_h,GAL3YJHH_p]=ttest2(data_output(2).values,data_output(1).values);
[GAL3BCHH_h,GAL3BCHH_h_p]=ttest2(data_output(3).values,data_output(1).values);

%Mann Whitney U-test (%Check syntax)
%[p,h,stats] = ranksum(data_output(2).values,data_output(1).values)
%[h,p,stats]=ranksum(data_output(3).values,data_output(1).values)

%%
strains={'RYC69';
    'RYD40';
    'RYB22'
    };
filename='SOK1_HH';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('../outputFigures/data_output_figure_SOK1HH','data_output');

[SOK1BC_h1,SOK1BC_p1]=ttest2(data_output(1).values,data_output(3).values);
[SOK1YJ_h2,SOK1YJ_p2]=ttest2(data_output(2).values,data_output(3).values);

%% CREATE LOG with all the numbers
%>>> EXPORT TO LOG
load('../outputFigures/log_results.mat');
cell2csv('../outputFigures/log_results.csv',log_results);











