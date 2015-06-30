function compute_ANOVA_to_compare_groups
%Determine if strains are statistically different using ANOVA for the
%comparison of groups

load('data_output_figure_1');
Strain_names={data_output.strain};

%YJM975,YJM978, YJM981
%QueryStrains_idx=[30:32];

%YOS128,YPS163,YPS606
%QueryStrains_idx=[34:36];

%Find indeces of matching strains
% ANOVA_Query_Strains={'DBVPG1373';
%     'DBVPG1106';
%     'L-1374';
%     'L-1528';
%     'BC187';
%     'YJM975';
%     'YJM978';
%     'YJM981';
%     'Bb32';
%     'M22';
%     'DBVPG1788';
%     'DBVPG6765';
%     'WE372'};

%Find indeces of matching strains
%Sake_strains={'UC5';'Y12-SGRP';'Y12-WashU';'Y9-WashU'};

%ANOVA_Query_Strains={'UC5';'Y12-SGRP';'Y12-WashU';'Y9-WashU'}

%% From haplotypes ORF only

ANOVA_Query_Strains={'YJM975','YJM978','YJM981'};
[p,tbl,stats]=compute_ANOVA(ANOVA_Query_Strains,data_output,Strain_names)


%ANOVA_Query_Strains={'BC187','DBVPG1373','DBVPG1788','DBVPG6765','L-1528'};

% ANOVA_Query_Strains={'DBVPG1106','L-1374','Bb32','UWOPS87-242.1','YS9'};
%ANOVA_Query_Strains={'DBVPG1106','L-1374','RM11_1A','UWOPS87_2421','YS9'};

%
%ANOVA_Query_Strains={'I-14','322134S','YS2'}
%ANOVA_Query_Strains={'I14','322134S','YS2'}

%
% ANOVA_Query_Strains={'273614N','YIIc17_E5'}
%
ANOVA_Query_Strains={'YPS163','T7','378604X','UWOPS83_787_3','YPS128','YPS606'}
[p,tbl,stats]=compute_ANOVA(ANOVA_Query_Strains,data_output,Strain_names)

%
% %SAME ORF but very different promoter?
%ANOVA_Query_Strains={'YJM421','REF','S288C','W303','YJM789'}
%ANOVA_Query_Strains={'YJM421','REF','S288c','W303','YJM789'}
%
% ANOVA_Query_Strains={'DBVPG6040','NCYC361'}
%
%ANOVA_Query_Strains={'UWOPS03_461_4','UWOPS05_217_3','UWOPS05_227_2'}

%ANOVA_Query_Strains={'Y9combined','Y9'}


%% From haplotypes PROMOTER + ORF

ANOVA_Query_Strains={'BC187','DBVPG1788','DBVPG6765','L-1528'};
[p,tbl,stats]=compute_ANOVA(ANOVA_Query_Strains,data_output,Strain_names)

ANOVA_Query_Strains={'YPS163','T7'};
[p,tbl,stats]=compute_ANOVA(ANOVA_Query_Strains,data_output,Strain_names)



