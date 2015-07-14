function subQueryStrainsMatrix = subSampleGeneticDistance(QueryStrains_matrix, matchingStrains)

% script to making a smaller phylogeny distance matrix

totalMatch=zeros(size(QueryStrains_matrix,1),2);
for i=1:length(matchingStrains)
    totalMatch(:,1)=totalMatch(:,1)+strcmpi(matchingStrains{i},QueryStrains_matrix(:,1));   
    totalMatch(:,2)=totalMatch(:,2)+strcmpi(matchingStrains{i},QueryStrains_matrix(:,2));
end

subQueryStrainsMatrix=QueryStrains_matrix(totalMatch(:,1)&totalMatch(:,2),:);

for i=1:length(subQueryStrainsMatrix)
    subQueryStrainsMatrix{i,5}=mean(mean(subQueryStrainsMatrix{i,4}));
end