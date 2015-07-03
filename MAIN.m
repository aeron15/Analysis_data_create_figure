function MAIN

%MAIN contains the highest level scripts

%Create an output folder if it does not exist
if ~exist('../outputFigures')
    mkdir('../outputFigures');
end

compute_setpoints_reference_BC187()
driver_main_figures()

%plot_compute_meiotic_segregants()
