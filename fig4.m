function fig4(diff_sp, names)

%% Figure 4. Natural Isolate ORF swaps into YJM978 

%strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*', 'RYD25*', 'RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*', 'RYD29*'};
strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYB92*', 'RYC72*', 'RYD25*', 'RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*'};

%Remvoed S288C
strains = {'RYC45*','RYC58*','RYC49*', 'RYC50*','RYC51*', 'RYC59_1*','RYC52*','RYC60*','RYC62*', 'RYC72*', 'RYD25*', 'RYD27*', 'RYD28*', 'RYD30*', 'RYD31*', 'RYB59*', 'RYB53*'};

%% PLOT NATURAL ISOLATE SETPOINTS AND ALLELE SWAP SETPOINTS

for iStrain = 1:length(strains)
        
    curr_strain = regexp(names, regexptranslate('wildcard', strains(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);

    sp_mean_swaps(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error_swaps(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
 
end

natural_isos = natural_isos_ref(strains);

for iStrain = 1:length(natural_isos)
    
    set_temp(:,iStrain) = strcmp(natural_isos{iStrain}, names);
    idx = find(set_temp(:,iStrain)==1);
    sp_mean(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
    
end

[num, loc] = sort(sp_mean_swaps);

k = figure('Position',[417   279   571   284]);

hold all;

errorbar(1:length(sp_mean),sp_mean(loc),std_error(loc),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);

errorbar(1:length(num),num,std_error_swaps(loc),'ok','MarkerFaceColor','red','MarkerSize',8);

lab = strain_name_conv(natural_isos);

xlim([0 length(lab)+1]);
xticklabel_rotate([1:length(lab)],45,lab(loc),'interpreter','tex');

Set_fig_RE(k,9,9,9);
ylim([-9 -3]);

filename = ['fig4_sorted_swaps.pdf'];
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

