function filtered_setpoints_normalized=extract_data(strains, setpoints_normalized)

names=setpoints_normalized(:,3);
counter=1;

for iStrain = 1:length(strains)
    
    curr_strain = regexp(names, regexptranslate('wildcard', strains(iStrain)));
    cs = cellfun(@isempty,curr_strain);
    
     idx = find(cs==0);
     
%      data_output(iStrain).strain=strains(iStrain);
%      data_output(iStrain).values=setpoints_normalized(idx);
     
     for i=1:length(idx)
        
         filtered_setpoints_normalized{counter,3}=strains(iStrain);
         
         filtered_setpoints_normalized{counter,1}=setpoints_normalized{idx(i),1};
         filtered_setpoints_normalized{counter,2}=setpoints_normalized{idx(i),2};
         
         counter=counter+1;
         
     end
%     sp_mean(iStrain) = mean(diff_sp(idx));
%     sp_std(iStrain) = std(diff_sp(idx));
%     std_error(iStrain) = sp_std(iStrain)/sqrt(length(idx)-1);
    
end

end