function [p,tbl,stats]=compute_ANOVA(ANOVA_Query_Strains,data_output,Strain_names)

counter=1;
for iStrain=1:length(ANOVA_Query_Strains)
    
    try
        QueryStrains_idx(counter)=find(strcmp(ANOVA_Query_Strains{iStrain},Strain_names));
        counter=counter+1;
    catch
        ANOVA_Query_Strains{iStrain}

    end
    
end


for iQueryStrain=1:length(QueryStrains_idx)
    
    
    QueryStrains_size(iQueryStrain)=length(data_output(QueryStrains_idx(iQueryStrain)).values);
end


ANOVA_Rows=max(QueryStrains_size);

ANOVA_Columns=length(QueryStrains_idx);

ANOVA_matrix=nan(ANOVA_Rows,ANOVA_Columns);

for iQueryStrain=1:length(QueryStrains_idx)
    
    Values=data_output(QueryStrains_idx(iQueryStrain)).values;
    
    for jValues=1:length(Values)
        
        ANOVA_matrix(jValues,iQueryStrain)=Values(jValues);
    end
    
    
end

[p,tbl,stats] =anova1(ANOVA_matrix);