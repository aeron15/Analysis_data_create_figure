function plot_correlation_highlight_points()

%PLOT_CORRELATION_HIGHLIGHT_POINTS plots specific data points

load('NaturalIsolates_GeneticDistance_CorrelationData')

%% Highlight specific strains

%Highlighted_Strains={'YPS128','YPS163','YPS606'};
%Highlighted_Strains={'YJM978','YJM975','YJM981'};
%Highlighted_Strains={'BC187','YJM978'};

%Highlighted_Strains={'YJM421','S288C'};

%Strains with maximal genetic distance
Highlighted_Strains={'IL-01','UWOPS03-461.4'};

Highlighted_Strains={'YJM978','DBVPG1373'}


[Highlight_x,Highlight_y]=highlight_in_correlation_plot(QueryStrains_matrix,QueryStrains_MeanGeneticDistance,Highlighted_Strains);


%% Plot correlation of genetic distance and the pairwise difference between set points

x=[QueryStrains_matrix{:,3}];

hfig=figure();
hold all;
plot(x,QueryStrains_MeanGeneticDistance,'.')
plot(Highlight_x,Highlight_y,'r.','MarkerSize',40)
xlabel('Log2 distance between set points')
ylabel('Genetic distance')

Set_fig_RE(hfig,14,14,14)

%% Find strains with maximal genetic distance ('IL-01'    'UWOPS03-461.4')

%find(QueryStrains_MeanGeneticDistance==max(QueryStrains_MeanGeneticDistance))
%QueryStrains_matrix(536,:)

%Find strain with maximal difference minimal genetic distance
%find(abs(x-3.639)<0.001)

%%
Strain_filename=[Highlighted_Strains{:}];
filename=['Correlation_setpoint_genetic_distance_highlight_strains' Strain_filename] ;
export_fig(filename, '-pdf','-transparent','-nocrop');
end
