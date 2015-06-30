function [AlleleReplacementBackgrounds_mean,AlleleReplacementBackgrounds_std]=compute_fold_difference_across_backgrounds(data_output,StrainsWithBC187Allele_names,StrainsWithYJM978Allele_names)

%DETERMINE_FOLD_DIFFERENCE_ACROSS_BACKGROUNDS 
%This assumes that strains are sorted equally in both Strains BC187 and
%YJM978
DataStrains_names={data_output.strain};

for iStrain=1:length(StrainsWithBC187Allele_names)
   
    %idx1=determine_index(DataStrains_names,StrainsWithBC187Allele_names{iStrain});
    %idx2=determine_index(DataStrains_names,StrainsWithYJM978Allele_names{iStrain});
    
    strain1=StrainsWithBC187Allele_names{iStrain}; 
    strain2=StrainsWithYJM978Allele_names{iStrain};
    [~,~,FoldDifferenceMean,~]=compute_fold_difference(data_output,strain1,strain2);
    FoldDifference_vector(iStrain)=FoldDifferenceMean;

end

%Compute men fold differnce and standard deviation

AlleleReplacementBackgrounds_mean=mean(FoldDifference_vector);
AlleleReplacementBackgrounds_std=std(FoldDifference_vector);

hfig=figure();
boxplot(FoldDifference_vector)
title('Distributions of fold differences between BC187 and YJM978 GAL3 alleles across strains')
Set_fig_RE(hfig,9,9,20);

filename='Box_Plot_of_fold differences';
export_fig(filename, '-pdf','-transparent','-nocrop');

end