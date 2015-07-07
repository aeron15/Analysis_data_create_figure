function [Highlight_x,Highlight_y]=highlight_in_correlation_plot(QueryStrains_matrix,QueryStrains_MeanGeneticDistance,Highlighted_Strains)

%HIGHLIGHT_IN_CORRELATION_PLOT hightlights specific points in the data set

counter=1;

for iStrain=1:length(Highlighted_Strains)
    
    for jStrain=iStrain+1:length(Highlighted_Strains)
        
        Strain1=Highlighted_Strains{iStrain};
        Strain2=Highlighted_Strains{jStrain};
        
        Strain_names1={QueryStrains_matrix{:,1}};
        Strain_names2={QueryStrains_matrix{:,2}};
        
        idx1=strcmp(Strain_names1,Strain1);
        idx2=strcmp(Strain_names2,Strain2);
        
        PairStrains_idx=find(idx1&idx2);
        
        if isempty(PairStrains_idx)
            
            idx1=strcmp(Strain_names1,Strain2);
            idx2=strcmp(Strain_names2,Strain1);
            
            PairStrains_idx=find(idx1&idx2);
            
        end
        %%
        
        if PairStrains_idx
            try
                Highlight_x(counter)=QueryStrains_matrix{PairStrains_idx,3};
                Highlight_y(counter)=QueryStrains_MeanGeneticDistance(PairStrains_idx);
                counter=counter+1;
            catch
                display('klajsdla')
            end
        end
    end
    
end

end