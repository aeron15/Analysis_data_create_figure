function Number_of_Groups=T_test_walking(data_output, loc)
%T_test walking
% Take the first strain and then t-test with the next strain

idx_groups=1;
idx_failed_test=1;
QueryStrains_ttest_result=[];

for iStrain=1:length(loc)-1
    idx=loc(iStrain);
    
    sample1=data_output(loc(iStrain)).values;
    sample2=data_output(loc(iStrain+1)).values;
    
    h=ttest2(sample1,sample2);
    QueryStrains_ttest_result(iStrain)=h;
    Number_of_Groups=sum(QueryStrains_ttest_result)+1;%Number of dividers or 1s is  


end

Strain_Names={data_output.strain};
Strain_Names=Strain_Names(loc);
%Get divider strains
%DividerStrains={Strain_Names{logical(QueryStrains_ttest_result)}};
end