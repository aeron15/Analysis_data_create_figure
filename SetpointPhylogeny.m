function SetpointPhylogeny(QueryStrains_MeanGeneticDistance, QueryStrains_matrix, data_output)

% The function SETPOINTPHYLOGENY is designed to create a phylogenetic tree
% from pairwise genetic distances of a given strains

% Modifed by KL 20150710 (to be function within
% compute_correlation_genetic_distance_setpoint_induction)

load('strainSources.mat');

%% Remove NaN values

idx_to_remove = isnan(QueryStrains_MeanGeneticDistance);
QueryStrains_MeanGeneticDistance(idx_to_remove) = [];

StrainNames = QueryStrains_matrix(:,1:2);
StrainNames(idx_to_remove,:) = [];

StrainNames_all = vertcat(StrainNames(1,1), StrainNames(1:31,2));

%% Construct Phylogentic tree

% tree = seqlinkage(QueryStrains_MeanGeneticDistance, 'average', StrainNames_all);
tree = seqneighjoin(QueryStrains_MeanGeneticDistance, 'equivar', StrainNames_all); %neighbor-joining tree

%% Categorize strains based on ecological niches - color bar
 
source_categories = unique({strainSource.source});
myColor = cbrewer('qual','Set1',6);

%% Construct Bar Chart

k = plot(tree);
set(gcf, 'color', [1 1 1]) ;
set(gcf, 'Position', [410 279 500 600])
box off;
filename = 'Setpoint_phylogeneticTree';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

Names_Order = get(tree, 'LeafNames');

hfig = figure('Position',[410   279   300   620]);

for iStrain = 1:length(Names_Order)
    idx_to_find_data = strcmp({data_output.strain}, Names_Order{iStrain});
    idx_to_find_source = strcmp({strainSource.strain}, Names_Order{iStrain});
    pos = length(Names_Order)+1-iStrain;
    
    values = data_output(idx_to_find_data).values;
    tSource=strainSource(idx_to_find_source).source;
    barh(pos, mean(values),'FaceColor',myColor(find(strcmpi(source_categories,tSource)),:),'BarWidth',.75);
    
    hold all;
end

ylim([0 length(Names_Order)+1]);
set(gca, 'YTick',1:length(Names_Order));
set(gca, 'YTickLabel',flip(Names_Order));

Set_fig_RE(hfig,9,9,9)

filename = 'Setpoint_phylogeny_bar_chart';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');




