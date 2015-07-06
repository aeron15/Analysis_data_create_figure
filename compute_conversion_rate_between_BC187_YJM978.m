function [YJM978bg_PercentConversion,BC187bg_PercentConversion,Error_conversion_YJ,Error_conversion_BC]=compute_conversion_rate_between_BC187_YJM978()

%Computes conversion rate and the error of this conversion using error
%propagation formulas

load('../outputFigures/data_output_figure_3.mat')


Strain_GAL3BC_YJ='GAL3-BC (YJM978)';
Strain_GAL3YJ_YJ='GAL3-YJM (YJM978)';
Strain_GAL3BC_BC='GAL3-BC (BC187) ';
Strain_GAL3YJ_BC='GAL3-YJM (BC187) ';

%Determine indeces of the strains
StrainsData_names={data_output.strain};

GAL3BC_YJ_idx=determine_index(StrainsData_names,Strain_GAL3BC_YJ);
GAL3YJ_YJ_idx=determine_index(StrainsData_names,Strain_GAL3YJ_YJ);
GAL3BC_BC_idx=determine_index(StrainsData_names,Strain_GAL3BC_BC);
GAL3YJ_BC_idx=determine_index(StrainsData_names,Strain_GAL3YJ_BC);

%Mean values and std for each strain
GAL3YJ_YJ_mean=nanmean(data_output(GAL3YJ_YJ_idx).values);
GAL3YJ_BC_mean=nanmean(data_output(GAL3YJ_BC_idx).values);
GAL3BC_YJ_mean=nanmean(data_output(GAL3BC_YJ_idx).values);
GAL3BC_BC_mean=nanmean(data_output(GAL3BC_BC_idx).values);

GAL3YJ_YJ_std=nanstd(data_output(GAL3YJ_YJ_idx).values);
GAL3YJ_BC_std=nanstd(data_output(GAL3YJ_BC_idx).values);
GAL3BC_YJ_std=nanstd(data_output(GAL3BC_YJ_idx).values);
GAL3BC_BC_std=nanstd(data_output(GAL3BC_BC_idx).values);


% Compute difference between endognenous alle and background
InterStrain_distance=abs(nanmean(data_output(GAL3BC_BC_idx).values)-nanmean(data_output(GAL3YJ_YJ_idx).values));

std_1=nanstd(data_output(GAL3BC_BC_idx).values);
std_2=nanstd(data_output(GAL3YJ_YJ_idx).values);
InterStrain_error=sqrt(std_1.^2+std_2.^2);


%Difference allele swaps

YJM978background_differences=abs(nanmean(data_output(GAL3YJ_YJ_idx).values)-nanmean(data_output(GAL3BC_YJ_idx).values));

std_1=nanstd(data_output(GAL3YJ_YJ_idx).values);
std_2=nanstd(data_output(GAL3BC_YJ_idx).values);
YJM978background_differences_error=sqrt(std_1.^2+std_2.^2);

BC187background_differences=abs(nanmean(data_output(GAL3BC_BC_idx).values)-nanmean(data_output(GAL3YJ_BC_idx).values));

std_1=nanstd(data_output(GAL3BC_BC_idx).values);
std_2=nanstd(data_output(GAL3YJ_BC_idx).values);
BC187background_differences_error=sqrt(std_1.^2+std_2.^2);

%Compute percent conversion strains
YJM978bg_PercentConversion=YJM978background_differences./InterStrain_distance;
BC187bg_PercentConversion=BC187background_differences./InterStrain_distance;

%% Compute error 
Error_conversion_YJ=sqrt((YJM978background_differences_error./YJM978background_differences).^2+(InterStrain_error./InterStrain_distance).^2).*YJM978bg_PercentConversion;
Error_conversion_BC=sqrt((BC187background_differences_error./BC187background_differences).^2+(InterStrain_error./InterStrain_distance).^2).*BC187bg_PercentConversion;

end