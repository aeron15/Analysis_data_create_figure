%COMPUTE_SETPOINTS_REFERENCE_BC187_galactose

% This computes the renormalized data for galactose
% Compute confidence interval

path_data='/Users/RenanEscalante/Dropbox/Phenotypic_diversity/var_analysis_data/20150623_data/GAL/'

load([path_data 'setpoints_normalized.mat']);

%Compute the 95% confidence interval for estimate

BC187_vals_vector=cell2mat(setpoints_normalized(:,1));

BC187_mean=nanmean(BC187_vals_vector);

BC187_std=nanstd(BC187_vals_vector);

% Standard deviation and length of the vector
%SEM = std(x)/sqrt(length(x));  
BC187_len=sum(~isnan(BC187_vals_vector));

BC187_SEM = BC187_std/sqrt(BC187_len);  

figure;
hist(BC187_vals_vector,20)
vline(BC187_mean,'green')

%% Compute the 10% above the mean and 10% below the mean

higher_bound= 1.1* BC187_mean;
lower_bound=0.89 * BC187_mean;
vline(higher_bound)
vline(lower_bound)
title ('BC187 measurements across replicates')

idx_to_remove=find(~(higher_bound * BC187_mean & BC187_vals_vector < lower_bound));
setpoints_normalized(idx_to_remove,:)=[];

%Create variable equivalent to dif_sp.mat

all_strains_vals_vector=cell2mat(setpoints_normalized(:,2));
all_strains_names=setpoints_normalized(:,3);

%% Get all the strains in the study

all_strains  = {'Y55*'; 'NCYC110*'; 'L_1528*'; 'DBVPG6044*';
    'Y12_SGRP*'; 'W303*'; 'i378604X*'; 'DBVPG1373*';
    'YIIc17_E5*'; 
    'CLIB324*'; 'NC_02*'; 'PW5*'; 'YS4*'; 
    'Y9_WashU*'; 'IL_01*';
    'YPS128*'; 'DBVPG1788*'; 
    'DBVPG1853*'; 'L_1374*'; 'DBVPG1106*'; 'YJM421*';
    'Bb32*'; 
    'YJM653*'; 'YPS1009*'; 
    'YJM975*';
    'FL100*'; 'i273614N*';
    'BC187*'; 'YJM978*';
    'RYC45*';'RYC58*';'RYC49*'; 'RYC50*';'RYC51*'; 'RYC59_1*';'RYC52*';'RYC60*';'RYC62*'; 'RYB92*'; 'RYC72*'}
    
%     'RY16*'; 'RYB53*'; 'RYB59*'; 'RYB65*'; 'RYB66*'; 'RYB28*';
%     'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';
%     'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*';
%     'RYC45*';'RYC58*';'RYC49*'; 'RYC50*';'RYC51*'; 'RYC59_1*';'RYC52*';'RYC60*';'RYC62*'; 'RYB92*'; 'RYC72*'; 
%     'RYD25*'; 'RYD27*'; 'RYD28*'; 'RYD30*'; 'RYD31*'; 'RYB59*'; 'RYB53*'};


rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);
filename='Fig_1_all_strains_galactose_titration';

[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

%% Get all the set points of induction of all the strains

%Removed CLIB215
all_strains  = {'Y55*'; 'NCYC110*'; 'L_1528*'; 'DBVPG6044*';
    'Y12_SGRP*'; 'W303*'; 'i378604X*'; 'DBVPG1373*';
    'YIIc17_E5*'; 
    'CLIB324*'; 'NC_02*'; 'PW5*'; 'YS4*'; 
    'Y9_WashU*'; 'IL_01*';
    'YPS128*'; 'DBVPG1788*'; 
    'DBVPG1853*'; 'L_1374*'; 'DBVPG1106*'; 'YJM421*';
    'Bb32*'; 
    'YJM653*'; 'YPS1009*'; 
    'YJM975*';
    'FL100*'; 'i273614N*';
    'BC187*'; 'YJM978*'};

rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);
filename='Fig_1_natural_isolates_galactose_titration';

[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);











    


