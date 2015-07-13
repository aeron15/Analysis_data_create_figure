function determine_GeneticDistance()

% The function DETERMINE_GENETICDISTANCE is designed to create a phylogenetic tree
% from pairwise genetic distances of a given strains

% Created by KL 20150710

%% Determine strains from the Cromie distance matrix

QueryStrains_counter=determine_strains_in_Cromie;

QueryStrains_genetic_distance=get_genetic_distance(Strain_1,Strain_2)


load('data_output_figure_1.mat');
load('NaturalIsolates_GeneticDistance_CorrelationData.mat');

%% Remove NaN values

idx_to_remove = isnan(QueryStrains_MeanGeneticDistance);
QueryStrains_MeanGeneticDistance(idx_to_remove) = [];

StrainNames = QueryStrains_matrix(:,1:2);
StrainNames(idx_to_remove,:) = [];

StrainNames_all = vertcat(StrainNames(1,1), StrainNames(1:34,2));

%% Construct Phylogentic tree

% tree = seqlinkage(QueryStrains_MeanGeneticDistance, 'average', StrainNames_all);
tree = seqneighjoin(QueryStrains_MeanGeneticDistance, 'equivar', StrainNames_all); %neighbor-joining tree

%% Categorize strains based on ecological niches - color bar
 
source_categories = unique({data_output.source});
myColor(1,:)=[1,0,0];
myColor(2,:)=[0,1,0];
myColor(3,:)=[0,0,1];
myColor(4,:)=[1,1,0];
myColor(5,:)=[1,0,1];

%% Construct Bar Chart

k = plot(tree);
% filename = 'Setpoint_phylogeny';
% export_fig(filename, '-pdf','-transparent','-nocrop');

Names_Order = get(tree, 'LeafNames');

k = figure;

for iStrain = 1:length(Names_Order)
    idx_to_find = strcmp({data_output.strain}, Names_Order{iStrain});
    pos = length(Names_Order)+1-iStrain;
    
    values = data_output(idx_to_find).values;
    tSource=data_output(idx_to_find).source;
    barh(pos, mean(values),'FaceColor',myColor(find(strcmpi(source_categories,tSource)),:),'BarWidth',.75);
    
    hold all;
end

ylim([0 length(Names_Order)+1]);
set(gca, 'YTick',1:35);
set(gca, 'YTickLabel',flip(Names_Order));

filename = 'Setpoint_bar_chart';
export_fig(filename, '-pdf','-transparent','-nocrop');




