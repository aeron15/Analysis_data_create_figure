function compute_correlation_galactose_glucose

%% Plot correlation of glucose and galactose

% load('data_output_figure_galactose_titration.mat')
% data_output_galactose=data_output;
% 
% load('data_output_figure_glucose_titration.mat')
% data_output_glucose=data_output;

% Plot the natural isolates

load('data_output_natural_isolates_glucose_titration')
data_output_glucose=data_output;

load('data_output_natural_isolates_galactose_titration')
data_output_galactose=data_output;

%Plot the data output

strains_glucose={data_output_glucose.strain};
strains_galactose={data_output_galactose.strain};

common_strains=intersect(strains_galactose,strains_glucose);
common_strains(15)=[];
%% Plot the correlation galactose and glucose for 20 strains
figure;
hold all;
xlim([-6.5 -2.5])
ylim([-9 -3])

for iStrain=1:length(common_strains)
    
    common_strains{iStrain}
    
    idx_galactose=find(strcmp(common_strains{iStrain}, strains_galactose));
    
    idx_glucose=find(strcmp(common_strains{iStrain}, strains_glucose));
    
    x1=nanmean(data_output_galactose(idx_galactose).values);
    x2=nanmean(data_output_glucose(idx_glucose).values);
    
    plot(x1,x2,'.')
    
    vector_x1(iStrain)=x1;
    vector_x2(iStrain)=x2;
    
end

xlabel('Galactose set point')
ylabel('Glucose set point')
axis square

% Check which strains have set points of induction that are off, check
% where the data for the supplementary material came from)

%IL-01
%Y12-SGRP
%YJM421

display('done')

