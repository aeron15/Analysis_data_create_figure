function MAIN

%MAIN contains the highest level scripts

%Create an output folder if it does not exist
if ~exist('../outputFigures')
    mkdir('../outputFigures');
end

%Analysis of data
% compute_setpoints_reference_BC187()
% 
% %Comparison of data in figure 1
% compute_ANOVA_to_compare_groups
% 
% driver_main_figures()

plot_compute_meiotic_segregants()

