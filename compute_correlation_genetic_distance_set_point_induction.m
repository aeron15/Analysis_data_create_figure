function [NaturalIsolates_correlation,x,QueryStrains_MeanGeneticDistance]=compute_correlation_genetic_distance_set_point_induction(data_output,loc)

%Computes the genetic distance and the set point of induction correlation
%between pairs of strains
%Sort the data
%Highlight specific strains on the plot
%Save strain names and get their location

data_output=data_output(loc);
strains={data_output.strain};

counter=1;

for iStrain=1:length(strains)
    
    for jStrain=iStrain+1:length(strains)%Makes fewer computations

        Strain1_mean_setpoint=nanmean(data_output(iStrain).values);
        Strain2_mean_setpoint=nanmean(data_output(jStrain).values);
        
        QueryStrains_Difference_between_setpoint=abs(Strain1_mean_setpoint-Strain2_mean_setpoint);
        
        QueryStrains_genetic_distance=get_genetic_distance(strains{iStrain},strains{jStrain});
        
        QueryStrains_matrix{counter,1}=strains{iStrain};
        QueryStrains_matrix{counter,2}=strains{jStrain};
        QueryStrains_matrix{counter,3}=QueryStrains_Difference_between_setpoint;
        QueryStrains_matrix{counter,4}=QueryStrains_genetic_distance;
                
        counter=counter+1;
    end
    
end
display('done')
%Plot correlation
%Average the distance for the cases where there is more than one distance
%computed
QueryStrains_MeanGeneticDistance=cellfun(@(x) nanmean(x(:)), QueryStrains_matrix(:,4), 'UniformOutput', false);
QueryStrains_MeanGeneticDistance=cell2mat(QueryStrains_MeanGeneticDistance);

%% Plot correlation of genetic distance and the pairwise difference between set points

x=[QueryStrains_matrix{:,3}];

hfig=figure();
plot(x,QueryStrains_MeanGeneticDistance,'.')
xlabel('Log2 distance between set points')
ylabel('Genetic distance')

Set_fig_RE(hfig,14,14,14)

%%
filename='NaturalIsolates_correlation_genetic distance_setpoint';
export_fig(filename, '-pdf','-transparent','-nocrop');

%%

NaturalIsolates_correlation=nancorr(x,QueryStrains_MeanGeneticDistance);

[x1, y]=remove_nan_rows(x,QueryStrains_MeanGeneticDistance');

[h,p]=corrcoef(x1,y);

save('NaturalIsolates_GeneticDistance_CorrelationData','QueryStrains_matrix','QueryStrains_MeanGeneticDistance')


end