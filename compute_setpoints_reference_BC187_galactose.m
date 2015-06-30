function compute_setpoints_reference_BC187_galactose()
%COMPUTE_SETPOINTS_REFERENCE_BC187_GALACTOSE

%% Load data

path_data='/Users/RenanEscalante/Dropbox/Phenotypic_diversity/var_analysis_data/20150623_data/GAL/';

load([path_data 'setpoints_normalized.mat']);

%% Filter data using BC187 set points. 

setpoints_normalized = filter_SetPointsNormalized(setpoints_normalized);

%% Create variable equivalent to dif_sp.mat

all_strains_vals_vector=cell2mat(setpoints_normalized(:,2));
all_strains_names=setpoints_normalized(:,3);

%Compute BC187 
BC187_vals_vector=cell2mat(setpoints_normalized(:,1));

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
    'RYC45*';'RYC58*';'RYC49*'; 'RYC50*';'RYC51*'; 'RYC59_1*';'RYC52*';'RYC60*';'RYC62*'; 'RYB92*'; 'RYC72*'};
    
%     'RY16*'; 'RYB53*'; 'RYB59*'; 'RYB65*'; 'RYB66*'; 'RYB28*';
%     'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';
%     'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*';
%     'RYC45*';'RYC58*';'RYC49*'; 'RYC50*';'RYC51*'; 'RYC59_1*';'RYC52*';'RYC60*';'RYC62*'; 'RYB92*'; 'RYC72*'; 
%     'RYD25*'; 'RYD27*'; 'RYD28*'; 'RYD30*'; 'RYD31*'; 'RYB59*'; 'RYB53*'};


rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);
filename='Fig_1_all_strains_galactose_titration';

[data_output,loc]=make_dot_plot(strains, all_strains_vals_vector, all_strains_names, filename);

data_output(1).values=nanmean(BC187_vals_vector);
data_output_galactose=data_output;
save('data_output_figure_galactose_titration','data_output_galactose');

filename = ['NaturalIsolatesStrains_galactose_setpoints.pdf'];
export_fig(filename, '-pdf','-transparent','-nocrop');


%% Get all the set points of induction of natural isolates

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

data_output_galactose=data_output;
save('data_output_natural_isolates_galactose_titration','data_output_galactose');




filename = ['AllStrain_galactose_setpoints.pdf'];
export_fig(filename, '-pdf','-transparent','-nocrop');














    


