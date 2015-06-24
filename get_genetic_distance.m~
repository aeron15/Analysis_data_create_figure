function QueryStrains_genetic_distance=get_genetic_distance(Strain_1,Strain_2)

%Determines the genetic distance between strains

load('Cromie_distances_2.mat')

%Check if replacing is needed

Strain_1=replace_strain(Strain_1);
Strain_2=replace_strain(Strain_2);

%Find index of the strain

curr_strain = regexp(names_strains, Strain_1);
cs = cellfun(@isempty,curr_strain);
Strain1_idx = find(cs==0);


curr_strain = regexp(names_strains, Strain_2);
cs = cellfun(@isempty,curr_strain);
Strain2_idx = find(cs==0);

%If the indeces exist for both strains
if ~(isempty(Strain1_idx) & isempty(Strain2_idx))
    
    QueryStrains_genetic_distance=patristic_distance_matrix(Strain1_idx,Strain2_idx);
else
    QueryStrains_genetic_distance=[];
end

end

function Strain=replace_strain(Strain)

if strcmp(Strain,'Bb32')
    
    Strain='RM11';
    
    
end

end

