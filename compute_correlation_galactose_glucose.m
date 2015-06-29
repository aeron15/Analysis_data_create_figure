function StrainsGalactoseGlucosenan_correlation=compute_correlation_galactose_glucose

%COMPUTE_CORRELATION_GALACTOSE_GLUCOSE computes the correlation of the mean
%set points of galactose and glucose

% Plot the natural isolates

load('data_output_natural_isolates_glucose_titration')
data_output_glucose=data_output;

load('data_output_natural_isolates_galactose_titration')


%Plot the data output

strains_glucose={data_output_glucose.strain};
strains_galactose={data_output_galactose.strain};

common_strains=intersect(strains_galactose,strains_glucose);
common_strains(15)=[];
%% Plot the correlation galactose and glucose for 20 strains

for iStrain=1:length(common_strains)
    
    %common_strains{iStrain};
    
    idx_galactose=find(strcmp(common_strains{iStrain}, strains_galactose));
    
    idx_glucose=find(strcmp(common_strains{iStrain}, strains_glucose));
    
    x1=nanmean(data_output_galactose(idx_galactose).values);
    x2=nanmean(data_output_glucose(idx_glucose).values);
    
    
    
    QueryStrain_SetPointGalactose(iStrain)=x1;
    QueryStrain_SetPointGlucose(iStrain)=x2;
    
    vector_x1(iStrain)=x1;
    vector_x2(iStrain)=x2;
    
end
%%

hfig=figure;
hold all;

plot(QueryStrain_SetPointGalactose,QueryStrain_SetPointGlucose,'.','MarkerSize',15);

% Plot the fitted line
[fittedX, fittedY]=compute_fit(QueryStrain_SetPointGalactose,QueryStrain_SetPointGlucose)
plot(fittedX, fittedY, 'r-', 'LineWidth', 3,'MarkerSize',15)


xlabel('Galactose set point')
ylabel('Glucose set point')

xlim([-6.5 -2.5])
ylim([-9 -3])

Set_fig_RE(hfig,9,9,9);
axis square

filename = ['Correlation_glucose_galactose_setpoints.pdf'];
export_fig(filename, '-pdf','-transparent','-nocrop');


%% Check which strains have set points of induction that are off, check
% where the data for the supplementary material came from)

%IL-01
%Y12-SGRP
%YJM421

%% Compute correlation coefficient

StrainsGalactoseGlucosenan_correlation=nancorr(QueryStrain_SetPointGalactose,QueryStrain_SetPointGlucose);
