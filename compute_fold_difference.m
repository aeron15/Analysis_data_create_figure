function [FoldDifferenceLowerBound,FoldDifferenceHigherBound,FoldDifferenceMean,ErrorFoldDifference]=compute_fold_difference(data_output,Strain1,Strain2)

%COMPUTE_FOLD_DIFFERENCE computes a higher and lower bound for the values
%of the set point of induction
StrainsData_names={data_output.strain};

Strain1_idx=determine_index(StrainsData_names,Strain1);
Strain2_idx=determine_index(StrainsData_names,Strain2);

Strain1_vals=data_output(Strain1_idx).values;
Strain2_vals=data_output(Strain2_idx).values;

Strain1_mean=nanmean(data_output(Strain1_idx).values);
Strain2_mean=nanmean(data_output(Strain2_idx).values);

Strain1_SEM=compute_standard_error(Strain1_vals);
Strain2_SEM=compute_standard_error(Strain2_vals);


Strain1_std=nanstd(Strain1_vals);
Strain2_std=nanstd(Strain2_vals);

%Compute lower bound

low1=Strain1_mean-Strain1_SEM;
low2=Strain2_mean-Strain2_SEM;
FoldDifferenceLowerBound=2.^(abs(low1-low2));

%Compute higher bound

high1=Strain1_mean+Strain1_SEM;
high2=Strain2_mean+Strain2_SEM;
FoldDifferenceHigherBound=2.^(abs(high1-high2));


FoldDifferenceMean=2.^abs(Strain1_mean-Strain2_mean);

ErrorFoldDifference=sqrt((Strain1_std).^2+(Strain2_std).^2);
%^BC187_YJM978_error_difference_Strains=sqrt((BC187_std_fig1).^2+(YJM978_std_fig1).^2);

%>>>>>EXPORT TO LOG
add_entry_log(['Fold Range difference between ' Strain1 ' and ' Strain2], FoldDifferenceMean);
add_entry_log(['Error fold difference between ' Strain1 ' and ' Strain2], ErrorFoldDifference);


end
