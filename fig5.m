function fig6(diff_sp, names)
%% Figure 6a. BC187/YJM978 allele swap in Natural Isolate background

strains_BC = {'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*'};
strains_YJM = {'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'};

%%
for iStrain = 1:length(strains_YJM)
        
    curr_strain = regexp(names, regexptranslate('wildcard', strains_YJM(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);

    sp_mean_YJM(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error_YJM(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
 
end

for iStrain = 1:length(strains_BC)
        
    curr_strain = regexp(names, regexptranslate('wildcard', strains_BC(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);

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

    sp_mean(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
    
end
%%
k = figure('Position',[417   313   400   200]);

errorbar(1:length(bc_sp),bc_sp,std_error_BC(loc),'ok','MarkerFaceColor','blue','MarkerSize',8);

hold all;

errorbar(1:length(yjm_sp),yjm_sp,std_error_YJM(loc),'ok','MarkerFaceColor','red','MarkerSize',8);

lab = strain_name_conv(natural_isos(loc));

ylim([-9 -3]);

xlim([0 length(lab)+1]);
xticklabel_rotate([1:length(lab)],45,lab,'interpreter','tex');

Set_fig_RE(k,9,9,9);
box off;
filename = ['fig5_both_alleles.pdf'];
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

%% plot setpoint gradients with error bars

strains = [strains_BC strains_BC];

close all;

%% PLOTS THE NATURAL ISOLATE, BC187, and YJM978 ALLELE STACKED

[sps, pos] = sort(sp_mean);

k = figure('Position',[417   313   400   200]);

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
box off;

filename = ['fig5_both_alleles_stacked.pdf'];
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

