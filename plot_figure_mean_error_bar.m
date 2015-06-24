function plot_figure_mean_error_bar(data_output)

%PLOT_FIGURE_MEAN_ERROR_BAR computes the mean of the data and standard error to plot data

%% compute all set point values

for iStrain=1:length(data_output)
    
    mean_data(iStrain)=mean(data_output(iStrain).values);
    standard_error(iStrain)=std(data_output(iStrain).values)./sqrt(length(data_output(iStrain).values));
    
end

%% Plot data

[meanDataSorted,idx]=sort(mean_data);

figure;
hold all;

for iStrain=1:length(data_output)
    
    plot(iStrain,meanDataSorted(iStrain),'o')
    %standard_error(iStrain)
    
end

%errorbar(meanDataSorted,standard_error(idx))
medianDataSort=median(meanDataSorted);
perc_75=prctile(meanDataSorted,75);
perc_25=prctile(meanDataSorted,25);

[counts,centers]=hist(meanDataSorted);

hfig=figure('Position',[ 328   324   870   271]);
hold all;
subplot(1,2,1)
errorbar(meanDataSorted,standard_error(idx))
hline(medianDataSort,'k')
hline(perc_75)
hline(perc_25)
subplot(1,2,2)
barh(centers,counts)
xlim([0 15])

namesStrains={data_output.strain};

hfig=figure('Position',[ 328   198   865   397]);
hold all;
positionVector1 = [0.1, 0.1, 0.6, 0.8];
subplot('Position',positionVector1)
errorbar(meanDataSorted,standard_error(idx),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8)
ylim([-5 2])
hline(medianDataSort,'k')
hline(perc_75)
hline(perc_25)
xticklabel_rotate(1:length(data_output),45,namesStrains(idx))
%set(gca,'box','off')

width=1;
positionVector2 = [0.64, 0.25, 0.1, 0.65];
subplot('Position',positionVector2)
p1=barh(centers,counts,width)
set(p1,'FaceColor',[0.5 0.5 0.5]);
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlim([0 15])
ylim([-5 2])

Set_fig_RE(hfig,9,9,20)

filename='Fig1_setpoints';