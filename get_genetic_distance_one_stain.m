function QueryStrains_genetic_distance=get_genetic_distance_one_stain(Strain_1)

%GET_GENETIC_DISTANCE Determines the genetic distance between strains
%It calls the function REPLACE_STRAIN

load('Cromie_distances_2.mat')

%Check if replacing is needed

Strain_1=replace_strain(Strain_1);

%Find index of the strain

curr_strain = regexp(CromieStrains_names, Strain_1);
cs = cellfun(@isempty,curr_strain);
Strain1_idx = find(cs==0);


Strain_1
curr_strain

end



