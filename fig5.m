function fig5(diff_sp, names)
%function fig5(diff_sp, names, glu_line)

% figure 5 of the GAL3 paper will contain data from the natural isolate
% allele swaps into YJM978

% created by KL 20150227

%% Figure 5. Natural Isolate ORF swaps into YJM978 

%strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*', 'RYD25*', 'RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*', 'RYD29*'};
strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*', 'RYD25*', 'RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*'};

% filename=['Paper_figs/fig5_natural_isos_YJM_20150403.pdf'];
% [tmp]=make_dot_plot_WC_renan_20150403(strains, diff_sp, names, filename);
% save('output_make_dotplot_fig4_swaps_in_YJ','tmp')

% filename=['Paper_figs/fig5_natural_isos_YJM_20150413.pdf'];
% fig_data = make_dot_plot_WC(strains, diff_sp, names, filename);

% filename=['Paper_figs/fig5_setpoint_corr'];
% fig5_corr(diff_sp, strains, names, filename);

% filename=['data_analysis/fig5_natural_isos_YJM.pdf'];
% make_dot_reps_WC(strains, diff_sp, names, filename);

% filename=['data_analysis/fig5_natural_isos_reps_NI_20150413.pdf'];
% make_dot_plot_sort_NI(strains, diff_sp, names, filename);

% close all;

%% plot setpoint gradients with error bars

% file = 'Paper_figs/S5_GLU_line_plot_';
% load_glu_data(glu_line, strains, file, 1);

% file = 'data_analysis/S5_GLU_line_plot_20150413';
% plot_glu_reps(glu_line, strains, file, 1);
% 
% close all;

%% PLOT NATURAL ISOLATE SETPOINTS AND ALLELE SWAP SETPOINTS

for iStrain = 1:length(strains)
        
    curr_strain = regexp(names, regexptranslate('wildcard', strains(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);

    hold all;
    sp_mean_swaps(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error_swaps(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
 
end

natural_isos = natural_isos_ref(strains);

for iStrain = 1:length(natural_isos)
    
    set_temp(:,iStrain) = strcmp(natural_isos{iStrain}, names);
    idx = find(set_temp(:,iStrain)==1);

    hold all;

    sp_mean(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
    
end

[num, loc] = sort(sp_mean_swaps);

k = figure('Position',[1   441   720   364]);

hold all;

errorbar(1:length(sp_mean),sp_mean(loc),std_error(loc),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);

errorbar(1:length(num),num,std_error_swaps(loc),'ok','MarkerFaceColor','red','MarkerSize',8);

lab = strain_name_conv(natural_isos);


xlim([0 length(lab)+1]);
xticklabel_rotate([1:length(lab)],45,lab(loc),'interpreter','tex');

Set_fig_RE(k,9,9,9);
ylim([-9 -3]);

%filename = ['Paper_figs/fig5_sorted_swaps.pdf'];
filename = ['fig5_sorted_swaps.pdf'];
export_fig(filename, '-pdf','-transparent','-nocrop');

%% plot histogram series
% 
% hists = load_hists_data(strains, glu_line, names);
% data_plot_1 = load_data(strains, glu_line, names, 2);
% 
% file = 'data_analysis/S5_hist_series_';
% plot_hist_series(hists, data_plot_1, file);
% 
% close all;