function NaturalIsolates_correlation=compute_correlation_genetic_distance_set_point_induction(data_output)

%Computes the genetic distance and the set point of induction correlation
%between pairs of strains

strains={data_output.strain};
counter=1;

for iStrain=1:length(strains)
    
    for jStrain=iStrain+1:length(strains)
        
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

%Plot correlation
%Average the distance for the cases where there is more than one distance
%computed

QueryStrains_mean_genetic_distance=cellfun(@(x) nanmean(x(:)), QueryStrains_matrix(:,4), 'UniformOutput', false);

QueryStrains_mean_genetic_distance=cell2mat(QueryStrains_mean_genetic_distance);
%% Plot correlation

x=[QueryStrains_matrix{:,3}];

figure(25)
plot(x,QueryStrains_mean_genetic_distance,'.')
xlabel('Difference set points')
ylabel('Genetic distance')

NaturalIsolates_correlation=nancorr(x,QueryStrains_mean_genetic_distance)
end