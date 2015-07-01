function plot_correlation_highlight_points()

%PLOT_CORRELATION_HIGHLIGHT_POINTS plots specific data points

load('NaturalIsolates_GeneticDistance_CorrelationData')

%% Highlight specific strains

Highlighted_Strains={'YPS128','YPS163','YPS606'};
%Highlighted_Strains={'YJM978','YJM975','YJM981'};

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

%%
%filename='Correlation_setpoint_genetic_distance_highlight_strains';
%export_fig(filename, '-pdf','-transparent','-nocrop');
end
