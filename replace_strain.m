function Strain=replace_strain(Strain)

%REPLACE_STRAIN replaces the name of some strains by other strain name that
%is equivalent for the scripts GET_GENETIC_DISTANCE and
%DETERMINE_STRAINS_IN_CROMIE

if strcmp(Strain,'Bb32')
    
    Strain='RM11';
    
end

if strcmp(Strain,'I-14')
    
    Strain='I14';
    
end

if strcmp(Strain,'I-14')
    
    Strain='I14';
    
end

if strcmp(Strain,'Y12-WashU')
    
    Strain='Y12';
    
end


if strcmp(Strain,'Y9-WashU')
    
    Strain='Y9';
    
end


if strcmp(Strain,'UWOPS87-242.1')
    
    Strain='UWOPS87-2421';
    
end

end
