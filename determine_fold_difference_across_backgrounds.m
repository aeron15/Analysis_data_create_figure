function determine_fold_difference_across_backgrounds(data_output,StrainsWithBC187Allele_names,StrainsWithYJM978Allele_names)

%DETERMINE_FOLD_DIFFERENCE_ACROSS_BACKGROUNDS 
%This assumes that strains are sorted equally in both Strains BC187 and
%YJM978
DataStrains_names={data_output.strain};

for iStrain=1:length(StrainsWithBC187Allele_names)
   
    %idx1=determine_index(DataStrains_names,StrainsWithBC187Allele_names{iStrain});
    %idx2=determine_index(DataStrains_names,StrainsWithYJM978Allele_names{iStrain});
    
    strain1=StrainsWithBC187Allele_names{iStrain}; 
    strain2=StrainsWithBC187Allele_names{iStrain};
    [FoldDifferenceLowerBound,FoldDifferenceHigherBound,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,strain1,strain2)

end


end