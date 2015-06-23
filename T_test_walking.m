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
    Number_of_Groups=sum(QueryStrains_ttest_result)+1;

    %Use the ones as separators between groups
    
%     if h==0
%         %If the test fails lump together the strains
%         failed_test(idx_failed_test)=idx;
%         idx_failed_test=idx_failed_test+1;
%         
%     end
%     
%     if h==1 %if test rejects null hypothesis
%         
%         if exist('failed_test')
%             
%             all_strain_cluster{idx_groups}=failed_test;
%             
%         else
%             
%             all_strain_cluster{idx_groups}=idx;
%             
%         end
%         
%         idx_groups=idx_groups+1;
%         idx_failed_test=1;
%         
%     end
%     
%     
% end
% 
%  all_strain_cluster{idx_groups}=failed_test;

end