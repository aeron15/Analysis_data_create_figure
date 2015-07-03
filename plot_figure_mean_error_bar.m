function plot_figure_mean_error_bar(data_output,varargin)

%PLOT_FIGURE_MEAN_ERROR_BAR computes the mean of the data and standard error to plot data
%% Parse parameters
p = inputParser;
addRequired(p,'data_output',@isstruct);
addParamValue(p,'pathOut','./',@isstr);
addParamValue(p,'file_append',date,@isstr);
addParamValue(p,'figure_size','default',@isstr);

parse(p,data_output,varargin{:});

pathOut=p.Results.pathOut;
file_append=p.Results.file_append;
figure_size=p.Results.figure_size;

%% Compute all set point values

for iStrain=1:length(data_output)
    
    mean_data(iStrain)=mean(data_output(iStrain).values);
    standard_error(iStrain)=std(data_output(iStrain).values)./sqrt(length(data_output(iStrain).values));
    
end

%% Plot data

[meanDataSorted,idx]=sort(mean_data);
medianDataSort=median(meanDataSorted);

namesStrains={data_output.strain};

switch figure_size
    case 'default'
        hfig=figure('Position',[ 328   198   865   397]);
        errorbar(meanDataSorted,standard_error(idx),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8)
        
    case 'HH'
        hfig=figure('Position',[328   337   358   258]);
        errorbar(meanDataSorted,standard_error(idx),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',6)
    
    case 'AlleleReplacement'
        hfig=figure('Position',[328   337   358   258]);
        errorbar(meanDataSorted,standard_error(idx),'ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',6)
    
        
end

hold all;

ylim([-9 -3])
xticklabel_rotate(1:length(data_output),45,namesStrains(idx));

set(gca,'box','off')
Set_fig_RE(hfig,9,9,20);

filename=[pathOut file_append 'Mean_error_bar_' ];
export_fig_specific_path(filename,'-pdf',  '-transparent', '-nocrop')
