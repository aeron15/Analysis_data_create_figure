function fig6(diff_sp, names)
%function fig6(diff_sp, names, glu_line)

% figure 5 of the GAL3 paper will contain data from the natural isolate
% allele swaps into YJM978

% created by KL 20150227

%% Figure 6a. BC187/YJM978 allele swap in Natural Isolate background

% strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD05*'; 'RYD12*'; 'RYD14*'; 'RYD53*'; 'RYB53_2'; 'RYD56*'; 'RYB65'};
% strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD06*'; 'RYD13*'; 'RYD15*'; 'RYD55*'; 'RYB59'; 'RYD58*'; 'RYB66'};

% strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD07*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'}; %'RYD50'; 'RYD53*'; 'RYD05*'
% strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD08*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'}; %'RYD52' 'RYD55*'; 'RYD06*'

% strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'; 'RYD50*'; 'RYD53*'; ''};
% strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*'};

strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'};
strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'};

%%

for iStrain = 1:length(strains_YJM)
        
    curr_strain = regexp(names, regexptranslate('wildcard', strains_YJM(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);

    hold all;
    sp_mean_YJM(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error_YJM(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
 
end

for iStrain = 1:length(strains_BC)
        
    curr_strain = regexp(names, regexptranslate('wildcard', strains_BC(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);

    hold all;
    sp_mean_BC(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error_BC(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
 
end

[num, loc] = sort(sp_mean_YJM);

bc_sp = sp_mean_BC(loc);
yjm_sp = sp_mean_YJM(loc);
natural_isos = natural_isos_ref(strains_YJM);

for iStrain = 1:length(natural_isos)
    
    set_temp(:,iStrain) = strcmp(natural_isos{iStrain}, names);
    idx = find(set_temp(:,iStrain)==1);

    hold all;

    sp_mean(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
    
end

k = figure('Position',[1   441   720   364]);

errorbar(1:length(bc_sp),bc_sp,std_error_BC(loc),'ok','MarkerFaceColor','blue','MarkerSize',8);

hold all;

errorbar(1:length(yjm_sp),yjm_sp,std_error_YJM(loc),'ok','MarkerFaceColor','red','MarkerSize',8);

lab = strain_name_conv(natural_isos(loc));

ylim([-9 -3]);

xlim([0 length(lab)+1]);
xticklabel_rotate([1:length(lab)],45,lab,'interpreter','tex');

Set_fig_RE(k,9,9,9);

filename = ['fig6_both_alleles.pdf'];
export_fig(filename, '-pdf','-transparent','-nocrop');

%% plot setpoint gradients with error bars

strains = [strains_BC strains_BC];
% file = 'Paper_figs/S6_GLU_line_plot_';
% load_glu_data(glu_line, strains, file, 1);

% file = 'data_analysis/S6_GLU_line_plot_';
% plot_glu_reps(glu_line, strains, file, 1);

close all;

%% PLOTS THE NATURAL ISOLATE, BC187, and YJM978 ALLELE STACKED

[sps, pos] = sort(sp_mean);

k = figure('Position',[1   441   720   364]);

errorbar(1:length(sps),sps,std_error(pos),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
hold all;
errorbar(1:length(sp_mean_BC),sp_mean_BC(pos),std_error_BC(pos),'ok','MarkerFaceColor','blue','MarkerSize',8);
errorbar(1:length(sp_mean_YJM),sp_mean_YJM(pos),std_error_YJM(pos),'ok','MarkerFaceColor','red','MarkerSize',8);

natural_isos = natural_isos_ref(strains_YJM(pos));
lab = strain_name_conv(natural_isos);

ylim([-9 -3]);

xlim([0 length(lab)+1]);
xticklabel_rotate([1:length(lab)],45,lab,'interpreter','tex');

Set_fig_RE(k,9,9,9);

filename = ['fig6_both_alleles_stacked.pdf'];
export_fig(filename, '-pdf','-transparent','-nocrop');


%% CORRELATION PLOTS

% strains = [strains_BC strains_YJM];
% filename=['Paper_figs/fig6_corr'];
% fig6_corr(diff_sp, strains, names, filename);

% filename=['Paper_figs/fig6_corr_BC'];
% fig_corr(diff_sp, strains_BC, names, filename);

% filename=['Paper_figs/fig6_corr_YJM'];
% fig_corr(diff_sp, strains_YJM, names, filename);

%% PLOT includes refererence strains -- just YJM

% k = figure('Position',[1   441   720   364]);
% 
% errorbar(1:length(yjm_sp),yjm_sp,std_error_YJM(loc),'ok','MarkerFaceColor','red','MarkerSize',8);
% 
% hold all;
% 
% errorbar(1:length(natural_isos),sp_mean,std_error,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
% 
% lab = strain_name_conv(natural_isos);
% 
% ylim([-5 1]);
% % ylab = {'3.1*10^{-2}', '6.3*10^{-2}', '1.3*10^{-1}', '2.5*10^{-1}', '0.5', '1'};
% ylab = {'0.031', '0.063', '0.130', '0.250', '0.50', '1.0', '2.0'};
% set(gca, 'Ytick',-5:1,'YTickLabel',ylab);
% 
% xlim([0 length(lab)+1]);
% xticklabel_rotate([1:length(lab)],45,lab,'interpreter','tex');
% 
% Set_fig_RE(figure,9,9,9);
% 
% filename = ['Paper_figs/fig6_YJM_ref.pdf'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
% 
% %% PLOT includes refererence strains -- just BC
% 
% k = figure('Position',[1   441   720   364]);
% 
% errorbar(1:length(bc_sp),bc_sp,std_error_BC(loc),'ok','MarkerFaceColor','blue','MarkerSize',8);
% 
% hold all;
% 
% errorbar(1:length(natural_isos),sp_mean,std_error,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
% 
% lab = strain_name_conv(natural_isos);
% 
% ylim([-5 1]);
% % ylab = {'3.1*10^{-2}', '6.3*10^{-2}', '1.3*10^{-1}', '2.5*10^{-1}', '0.5', '1'};
% ylab = {'0.031', '0.063', '0.130', '0.250', '0.50', '1.0', '2.0'};
% set(gca, 'Ytick',-5:1,'YTickLabel',ylab);
% 
% xlim([0 length(lab)+1]);
% xticklabel_rotate([1:length(lab)],45,lab,'interpreter','tex');
% 
% Set_fig_RE(figure,9,9,9);
% 
% filename = ['Paper_figs/fig6_BC_ref.pdf'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
% 
% %% PLOT includes refererence strains -- ALL
% 
% k = figure('Position',[1   441   720   364]);
% 
% errorbar(1:length(yjm_sp),yjm_sp,std_error_YJM(loc),'ok','MarkerFaceColor','red','MarkerSize',8);
% 
% hold all;
% 
% errorbar(1:length(natural_isos),sp_mean,std_error,'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
% 
% errorbar(1:length(bc_sp),bc_sp,std_error_BC(loc),'ok','MarkerFaceColor','blue','MarkerSize',8);
% 
% lab = strain_name_conv(natural_isos);
% 
% ylim([-5 1]);
% % ylab = {'3.1*10^{-2}', '6.3*10^{-2}', '1.3*10^{-1}', '2.5*10^{-1}', '0.5', '1'};
% ylab = {'0.031', '0.063', '0.130', '0.250', '0.50', '1.0', '2.0'};
% set(gca, 'Ytick',-5:1,'YTickLabel',ylab);
% 
% xlim([0 length(lab)+1]);
% xticklabel_rotate([1:length(lab)],45,lab,'interpreter','tex');
% 
% Set_fig_RE(figure,9,9,9);
% 
% filename = ['Paper_figs/fig6_all.pdf'];
% export_fig(filename, '-pdf','-transparent','-nocrop');
% 
% 
% 
% % k = figure('Position',[1   441   720   364]);
% % 
% % filename=['Paper_figs/fig6_natural_isos_YJM_BC.pdf'];
% % make_dot_plot_WC(strains, diff_sp, names, filename);
% % 
% % filename=['data_analysis/fig6_natural_isos_YJM_BC.pdf'];
% % make_dot_reps_WC(strains, diff_sp, names, filename);
% % 
% % close all;
% % 
% % %% plot setpoint gradients with error bars
% % 
% % data_plot = load_data(strains, glu_line, names, 1);
% % 
% % file = 'Paper_figs/S6_GLU_line_plot_';
% % load_glu_data(data_plot, file);
% % 
% % file = 'data_analysis/S6_GLU_line_plot_';
% % plot_glu_reps(data_plot, file);
% % 
% % close all;
% % 
% % %% plot histogram series
% % 
% % hists = load_hists_data(strains, glu_line, names);
% % data_plot_1 = load_data(strains, glu_line, names, 2);
% % 
% % file = 'data_analysis/S6_hist_series_';
% % plot_hist_series(hists, data_plot_1, file);
% % 
% % close all;