function plot_figure_mean_error_bar_side_histogram(data_output,varargin)

%PLOT_FIGURE_MEAN_ERROR_BAR computes the mean of the data and standard error to plot data

%% Parse parameters
p = inputParser;
addRequired(p,'data_output',@isstruct);
addParamValue(p,'pathOut','./',@isstr);
addParamValue(p,'file_append',date,@isstr);

parse(p,data_output,varargin{:});

pathOut=p.Results.pathOut;
file_append=p.Results.file_append;

%% compute all set point values

for iStrain=1:length(data_output)
    
    mean_data(iStrain)=mean(data_output(iStrain).values);
    standard_error(iStrain)=std(data_output(iStrain).values)./sqrt(length(data_output(iStrain).values));
    
end

%% Plot data

[meanDataSorted,idx]=sort(mean_data);

medianDataSort=median(meanDataSorted);
perc75=prctile(meanDataSorted,75);
perc25=prctile(meanDataSorted,25);

%Determine number of strains contained between the 25th and 75th percentile

PercentBetweenPerc25Perc75=sum(perc25<mean_data & mean_data<perc75)./length(mean_data);

[counts,centers]=hist(meanDataSorted);

namesStrains={data_output.strain};

hfig=figure('Position',[ 328   198   865   397]);
hold all;
positionVector1 = [0.1, 0.1, 0.6, 0.8];
subplot('Position',positionVector1)
errorbar(meanDataSorted,standard_error(idx),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8)
ylim([-9 -3])
hline(medianDataSort,'k')
hline(perc75)
hline(perc25)
xticklabel_rotate(1:length(data_output),45,namesStrains(idx))
set(gca,'box','off')

width=1;
positionVector2 = [0.64, 0.25, 0.1, 0.65];
subplot('Position',positionVector2)
p1=barh(centers,counts,width)
set(p1,'FaceColor',[0.5 0.5 0.5]);
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlim([0 15])
ylim([-9 -3])

Set_fig_RE(hfig,9,9,20)

filename=[pathOut 'Mean_error_bar_side_histogram' file_append];
export_fig(filename,'-pdf',  '-transparent', '-nocrop')
