function [data_output,loc]=make_dot_plot(strains, diff_sp, names, filename)

%MAKE_DOT_PLOT  makes a dot plot with errorbars where the x-axis is
%the name of the strains and the y axis is the set point of induction

lab = strain_name_conv(strains);

for iStrain = 1:length(strains)
    
    curr_strain = regexp(names, regexptranslate('wildcard', strains(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);
    
    data_output(iStrain).strain=lab{iStrain};
    data_output(iStrain).values=(diff_sp(idx));
    
    sp_mean(iStrain) = mean(diff_sp(idx));
    sp_std(iStrain) = std(diff_sp(idx));
    std_error(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
    
end

%%

[~, loc] = sort(sp_mean);
k = figure('Position',[1   441   720   364]);

strains_new = strains(loc);

for iStrain = 1:length(strains_new)
    
    curr_strain = regexp(names, regexptranslate('wildcard', strains_new(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
    idx = find(cs==0);
    
    hold all;
    
    sp_mean = diff_sp(idx);
    
    if isempty(sp_mean)
        sp_mean=nan;
    end
    plot(iStrain,sp_mean,'ok','MarkerSize',10,'MarkerFaceColor',[0.5 0.5 0.5]);
    hold all;
    plot(iStrain,mean(sp_mean),'ok','MarkerSize',10,'MarkerFaceColor',rgb('DarkOrange'));
    
end

format_fig(k, strains, loc)
export_fig(filename, '-pdf','-transparent','-nocrop');

%% Save data into a file
%save(['output_' filename],'data_output','loc')

