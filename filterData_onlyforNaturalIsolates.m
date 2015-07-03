function setpoints_normalized=filterData_onlyforNaturalIsolates(setpoints_normalized)

%filterData_onlyforNaturalIsolates
% Keep only the data in the paper

all_strains  = {'Y55*'; 'NCYC110*'; 'L_1528*'; 'DBVPG6044*';
    'Y12_SGRP*'; 'W303*'; 'i378604X*'; 'DBVPG1373*';
    'YIIc17_E5*'; 'UWOPS87_2421*'; 'YPS163*'; 'CLIB215*';
    'CLIB324*'; 'NC_02*'; 'PW5*'; 'YS4*'; 
    'T7*'; 'Y9_WashU*'; 'UWOPS03_4614*'; 'IL_01*';
    'M22*'; 'DBVPG6765*'; 'YPS128*'; 'DBVPG1788*'; 
    'DBVPG1853*'; 'L_1374*'; 'DBVPG1106*'; 'YJM421*';
    'Bb32*'; 'YJM428*'; 'UWOPS05_2272*'; 'DBVPG6040*';
    'YJM653*'; 'UC5*'; 'YPS1009*'; 'CLIB382*';
    'WE372*'; 'YJM975*'; 'I_14*'; 'YJM981*';
    'Y12_WashU*'; 'FL100*'; 'i273614N*'; 'YPS606*';
    'BC187*'; 'YJM978*'; 'S288C*';'RY16*'; 'RYB53*'; 'RYB59*'; 'RYB65*'; 'RYB66*'; 'RYB28*';
    'RYD42*'; 'RYD01*'; 'RYD03*'; 'RYD12*'; 'RYD14*'; 'RYB65*'; 'RYB53*';
    'RYB89*'; 'RYD02*'; 'RYD04*'; 'RYD13*'; 'RYD15*'; 'RYB66*'; 'RYB59*'; 'RYD52*'; 'RYD55*'; 'RYD06*';
    'RYC45*';'RYC58*';'RYC49*'; 'RYC50*';'RYC51*'; 'RYC59_1*';'RYC52*';'RYC60*';'RYC62*'; 'RYB92*'; 'RYC72*'; 
    'RYD25*'; 'RYD27*'; 'RYD28*'; 'RYD30*'; 'RYD31*'; 'RYB59*'; 'RYB53*';'RYB22';'RYB23';'RYB24';'RYC69';'RYD40';'RYB22'};

rm_strains = {'YIIc17_E5*'; 'i273614N*'; 'i378604X*'; 'YS4*'; 'NCYC110*'; 'Y55*'; 'PW5*'; 'DBVPG6044*'; 'W303*'; 'UWOPS05_2272*'};

strains = setdiff(all_strains, rm_strains);


%% Keeo only the ones that match all_strains

setpoints_normalized=extract_data(strains, setpoints_normalized);

%%

BC187_vals_vector=cell2mat(setpoints_normalized(:,1));
BC187_vals_vector(isnan(BC187_vals_vector))=[];

BC187_mean=nanmean(BC187_vals_vector);
BC187_median=median(BC187_vals_vector);

BC187_std=nanstd(BC187_vals_vector);
BC187_len=sum(~isnan(BC187_vals_vector));
BC187_SEM = BC187_std/sqrt(BC187_len);  %SEM = std(x)/sqrt(length(x));  

%%

higher_bound=BC187_median-BC187_std.*2;
lower_bound=BC187_median+BC187_std.*2;
idx_to_remove=(BC187_vals_vector<higher_bound | BC187_vals_vector>lower_bound);

%BC187_vals_vector(idx_to_remove)=[];

% %Cases to be removed
% %setpoints_normalized(idx_to_remove,:)
% 
removedStrains_number=sum(idx_to_remove);
setpoints_normalized(idx_to_remove,:)=[];

%% Remove instances where there is only one replicate of the data
strains_names=setpoints_normalized(:,3);
strain_nameConversion=short_strain_name(strains_names)';%From names in the data
strainsPaper_names=unique(strain_nameConversion);%This are unique identifiers used for labels in the paper

%Identify unique
strainsOneReplicate_counter=1;

for iStrain=1:length(strain_nameConversion)
    
    replicate_number=sum(strcmp(strain_nameConversion{iStrain},strain_nameConversion));%identify the number of replicates for that strain
    
    if replicate_number ==1%save idx for removal
        StrainsOneReplicate_idx(strainsOneReplicate_counter)=iStrain;
        strainsOneReplicate_counter=strainsOneReplicate_counter+1;
    end
    
end



    

strainsOneReplicate_number=length(StrainsOneReplicate_idx);
%Remove strains with one replicate
setpoints_normalized(StrainsOneReplicate_idx,:)=[];

%%
plot_hist_BC187_vals(BC187_vals_vector,higher_bound,lower_bound)

end