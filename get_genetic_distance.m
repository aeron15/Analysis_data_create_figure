function QueryStrains_genetic_distance=get_genetic_distance(Strain_1,Strain_2)

%GET_GENETIC_DISTANCE Determines the genetic distance between strains
%It calls the function REPLACE_STRAIN

load('Cromie_distances_2.mat')

CromieStrains_names(124) = {'BC186_Group_2_1.00'}; %change name so not called with YPS606
CromieStrains_names(119) = {'T3r_Group_5_0.85'}; %change name so not called with T7
CromieStrains_names(173) = {'T3g_Group_5_0.51'}; %change name so not called with T7

%Check if replacing is needed

Strain_1=replace_strain(Strain_1);
Strain_2=replace_strain(Strain_2);

%Find index of the strain

curr_strain = regexp(CromieStrains_names, Strain_1);
cs = cellfun(@isempty,curr_strain);
Strain1_idx = find(cs==0);


curr_strain = regexp(CromieStrains_names, Strain_2);
cs = cellfun(@isempty,curr_strain);
Strain2_idx = find(cs==0);

%If the indeces exist for both strains
if ~(isempty(Strain1_idx) & isempty(Strain2_idx))
    
    QueryStrains_genetic_distance=patristic_distance_matrix(Strain1_idx,Strain2_idx);
else
    QueryStrains_genetic_distance=[];
end

end


