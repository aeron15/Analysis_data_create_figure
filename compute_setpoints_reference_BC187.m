function compute_setpoints_reference_BC187()

%COMPUTE_SETPOINTS_REFERENCE_BC187

%Input should be the set points normalized in data
%Compute confidence interval for the true value of the set point of BC187
%Only keep measurements that fall in that interval

%Replot figure 1, with replicates

%I would like to see a comparison between of the
%relative setpoint (difference between setpoints) and then get the actual
%setpoint (perhaps look at rankings), but this is not necessary

%This should be a function that takes the mean of the replicates per
%strain and classifies strains based on tahat
%Perform a statistical test of different categories,that is, perform t test
%in order when the t test fails to rejec the null hypothesis that defines a
%new group


%This type of normalization has the advantage that it uses a more natural y
%axis. Also, it might solve some of the issues with the variability across
%experiments on figure 3 ( the comparison of allele swaps).

%% Compute confidence interval


path_data='/Users/RenanEscalante/Dropbox/Phenotypic_diversity/var_analysis_data/20150623_data/GLU/'

%path_data='20150623_data/GLU/'

load([path_data 'setpoints_normalized.mat']);

%Compute the 95% confidence interval for estimate

BC187_vals_vector=cell2mat(setpoints_normalized(:,1));

BC187_mean=nanmean(BC187_vals_vector);

BC187_std=nanstd(BC187_vals_vector);

% Standard deviation and length of the vector
%SEM = std(x)/sqrt(length(x));  
BC187_len=sum(~isnan(BC187_vals_vector));

BC187_SEM = BC187_std/sqrt(BC187_len);  

%ts = tinv([0.025  0.975],BC187_len-1);      % T-Score

ts = tinv([0.01  0.99],BC187_len-1);      % T-Score

CI = BC187_mean + ts*BC187_SEM;

%%
figure;
hist(BC187_vals_vector,20)
vline(BC187_mean,'green')
%vline(CI(1))
%vline(CI(2))

% Filter data based on if the BC187 set point falls in the 95% confidence
% interval

%CI(1)<=BC187_vals_vector 
%CI(2)>=BC187_vals_vector

%% Compute the 10% above the mean and 10% below the mean

higher_bound= 1.1* BC187_mean;
lower_bound=0.89 * BC187_mean;
vline(higher_bound)
vline(lower_bound)
title ('BC187 measurements across replicates')

idx_to_remove=find(~(higher_bound * BC187_mean & BC187_vals_vector < lower_bound));
setpoints_normalized(idx_to_remove,:)=[];

%%
%Create variable equivalent to dif_sp.mat

all_strains_vals_vector=cell2mat(setpoints_normalized(:,2));
all_strains_names=setpoints_normalized(:,3);

%% Figure 1. Setpoints of all strains in study

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

average_perc_difference_replicates=compute_percent_difference_between_replicates(data_output);

save('data_output_figure_glucose_titration','data_output');


%% Figure 1. Setpoints of natural isolates

%Removed CLIB215
all_strains  = {'Y55*'; 'NCYC110*'; 'L_1528*'; 'DBVPG6044*';
    'Y12_SGRP*'; 'W303*'; 'i378604X*'; 'DBVPG1373*';
    'YIIc17_E5*'; 'UWOPS87_2421*'; 'YPS163*';
    'CLIB324*'; 'NC_02*'; 'PW5*'; 'YS4*'; 
    'T7*'; 'Y9_WashU*'; 'UWOPS03_4614*'; 'IL_01*';
    'M22*'; 'DBVPG6765*'; 'YPS128*'; 'DBVPG1788*'; 
    'DBVPG1853*'; 'L_1374*'; 'DBVPG1106*'; 'YJM421*';
    'Bb32*'; 'YJM428*'; 'UWOPS05_2272*'; 'DBVPG6040*';
    'YJM653*'; 'UC5*'; 'YPS1009*'; 'CLIB382*';
    'WE372*'; 'YJM975*'; 'I_14*'; 'YJM981*';
    'Y12_WashU*'; 'FL100*'; 'i273614N*'; 'YPS606*';
    'BC187*'; 'YJM978*'; 'S288C*'};

rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);
filename='Fig_1_natural_isolates_20150609';

%make_dot_plot_WC_renan_20150403(strains, all_strains_vals_vector, all_strains_names, filename);
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
save('data_output_figure_1','data_output');
save('data_output_natural_isolates_glucose_titration','data_output');

% Average standard deviation for paper
average_standard_deviation=compute_average_standard_deviation(data_output);

%Correlation genetic distance and set point of induction
NaturalIsolates_correlation=compute_correlation_genetic_distance_set_point_induction(data_output);

Number_of_Groups=T_test_walking(data_output, loc);

% Compute difference between BC187 and YJM978
BC187_idx=find(strcmp({data_output.strain},'BC187'));
YJM978_idx=find(strcmp({data_output.strain},'YJM978'));

2.^abs(mean(data_output(BC187_idx).values)-mean(data_output(YJM978_idx).values))

BC187_std_fig1=std(data_output(BC187_idx).values);
YJM978_std_fig1=std(data_output(YJM978_idx).values);

BC187_YJM978_error_difference_Strains=sqrt((BC187_std_fig1).^2+(YJM978_std_fig1).^2);
%Range of the set points on figure 1
YJM421_last_loc_range=2.^abs(mean(data_output(loc(1)).values)-mean(data_output(loc(end)).values));

%Range of the strains on figure 4
YJM421_DBPVG1373_range=2.^abs(mean(data_output(27).values)-mean(data_output(6).values));

%% Figure 3
%strains = {'RY16*', 'RYB53*', 'RYB59*', 'RYB65*', 'RYB66*', 'RYB28*'};
strains = {'RYB53*', 'RYB59*', 'RYB65*', 'RYB66*'};

filename='Fig_3_allele_swaps_20150609';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

%compute_conversion_rate

InterStrain_distance=abs(mean(data_output(2).values)-mean(data_output(3).values));

data_output(1).values=data_output(1).values([1:3 5]);

YJM978background_differences=abs(mean(data_output(2).values)-mean(data_output(1).values))

BC187background_differences=abs(mean(data_output(3).values)-mean(data_output(4).values));


YJM978background_differences./InterStrain_distance
BC187background_differences./InterStrain_distance

%% Figure 3 test heterologos locus effect
strains = {'RY16*', 'RYB53*', 'RYB59*', 'RYB65*', 'RYB66*', 'RYB28*'};

filename='Fig_3_allele_swaps';
[fig3,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

%compute_conversion_rate

InterStrain_distance=abs(mean(fig3(2).values)-mean(fig3(3).values));

fig3(2).values=fig3(2).values([1:3 5]);

YJM978background_differences=abs(mean(fig3(2).values)-mean(fig3(1).values))

BC187background_differences=abs(mean(fig3(3).values)-mean(fig3(4).values));


YJM978background_differences./InterStrain_distance
BC187background_differences./InterStrain_distance

BC_het_effect=mean([fig3(4).values;fig3(6).values]);
YJ_het_effect=mean([fig3(1).values;fig3(3).values]);

Delta_GAL3=BC_het_effect-YJ_het_effect;

set_1=[fig3(1).values;fig3(6).values-Delta_GAL3];%haploids
set_2=[fig3(3).values;fig3(4).values-Delta_GAL3];

[h,p,ci,stats]=ttest2(set_1,set_2);


%% Figure 4. Natural Isolate ORF swaps into YJM978 

strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*','RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*'};
filename='Fig_4_YJ_bg_Diff_alleles';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);
T_test_walking(data_output, loc);

%Range of variation
2.^abs(mean(data_output(loc(1)).values)-mean(data_output(loc(end)).values))

save('data_output_figure_4','data_output');

[Correlation_Coefficient,P_Value]=compute_correlation_natural_isolates_allele_replacements;


%% Figure 5 BC
%strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'; 'RYD50*'; 'RYD53*'};
strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'};

filename='Fig5_BC';
make_dot_plot(strains_BC, all_strains_vals_vector, all_strains_names, filename);
%% Figure 5

strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*'};
filename='Fig5_YJ_alelle';
make_dot_plot(strains_YJM, all_strains_vals_vector, all_strains_names, filename);
display('done')


%% Figure 5 BC
strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*'};

filename='Fig5_BC_YJ';
make_dot_plot(strains_BC, all_strains_vals_vector, all_strains_names, filename);

%% SI figures. Hemizygous hybrid strains

strains={'RYB22';
'RYB23';
'RYB24'}
%'RYB42'}
filename='GAL3_HH';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

[h,p]=ttest2(data_output(2).values,data_output(1).values)
[h,p]=ttest2(data_output(3).values,data_output(1).values)

%%
strains={'RYC69';
'RYD40';
'RYB22'
};
filename='SOK1_HH';
[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

[h,p]=ttest2(data_output(1).values,data_output(3).values)
[h,p]=ttest2(data_output(2).values,data_output(3).values)








    


