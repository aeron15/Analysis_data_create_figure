function generateDistPhenotypeWar(subQueryStrainsMatrix, natural_isolate_matrix, NaturalIsolates_correlation)

% Generate correlation of traits from Warringer Data based on strains
% included in our study

matchingStrains = natural_isolate_matrix(1,:);
aa=cellfun(@str2num,natural_isolate_matrix(4:end,:),'UniformOutput',false);
aaa=cellfun(@isempty,aa);
for iPos=find(aaa)'
    aa{iPos}=nan;
end
subNatIsoMatWarNum=aa;
%%
iCount=0;
pairTraits=[];
for iPair=1:size(subQueryStrainsMatrix,1)
    iStrain1=subQueryStrainsMatrix{iPair,1};
    iStrain2=subQueryStrainsMatrix{iPair,2};
    posStrain1=strcmpi(iStrain1,natural_isolate_matrix(1,:));
    posStrain2=strcmpi(iStrain2,natural_isolate_matrix(1,:));
    distVector=abs([subNatIsoMatWarNum{4:end,posStrain1}]-[subNatIsoMatWarNum{4:end,posStrain2}]);
    iCount=iCount+1;
    pairTraits(iCount,:)=distVector;
end
%%
iCount=0;
pairALLTraits=[];
for iPair=1:size(subQueryStrainsMatrix,1)
    iStrain1=subQueryStrainsMatrix{iPair,1};
    iStrain2=subQueryStrainsMatrix{iPair,2};
    posStrain1=strcmpi(iStrain1,natural_isolate_matrix(1,:));
    posStrain2=strcmpi(iStrain2,natural_isolate_matrix(1,:));
    tmp=[subNatIsoMatWarNum{4:end,posStrain1}];
    tmp2=[subNatIsoMatWarNum{4:end,posStrain2}];
    removeThese=[find(isnan(tmp)),find(isnan(tmp2))];
    tmp(removeThese)=[];
    tmp2(removeThese)=[];
    distVector=corrcoef(tmp,tmp2);
    distVector=distVector(1,2);
    iCount=iCount+1;
    pairALLTraits(iCount,:)=distVector;
end
%%
myCorr=[];
for iCondition=1:size(pairTraits,2)
    a=corrcoef(pairTraits(:,iCondition),[subQueryStrainsMatrix{:,5}]);
    myCorr(iCondition,1)=a(1,2);
end

%% Plot histogram of Warringer Correlation

hfig=figure();

hist(myCorr);
line([NaturalIsolates_correlation, NaturalIsolates_correlation], ylim);
xlabel('R-values for correlation of Warringer Traits');

Set_fig_RE(hfig,9,9,9);
filename='WarringerTraits_genetic_distance_histogram';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');

%% Plot boxplot
boxplot(myCorr);
refline(0, NaturalIsolates_correlation);
title('R-values for correlation of Warringer Traits');

filename='WarringerTraits_genetic_distance_boxplot';
export_fig_specific_path(filename, '-pdf','-transparent','-nocrop');


%% Make a list of the sources of each strain
iCount=0;
for iStrain=matchingStrains
    iCount=iCount+1;
    posStrain=find(strcmpi(iStrain,natural_isolate_matrix(1,:)));
    infoStrain(iCount,1:2)=natural_isolate_matrix(2:3,posStrain);
end

