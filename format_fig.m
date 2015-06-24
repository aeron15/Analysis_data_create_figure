function format_fig_KL(figure, strains, loc)

% The FORMAT_FIG_KL function formats setpoint plots for the GAL3 paper

% Inputs are the figure handle, location of data points, and the strain
% names

% Created by KL 20150302

%% Apply plot formating

lab = strain_name_conv(strains);

ylim([-9.5 -3]);

xlim([0 length(lab)+1]);
xticklabel_rotate([1:length(lab)],45,lab(loc),'interpreter','tex');

Set_fig_RE(figure,9,9,9);
