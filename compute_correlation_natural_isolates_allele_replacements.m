function [Correlation_Coefficient,P_Value]=compute_correlation_natural_isolates_allele_replacements()

load('data_output_figure_1.mat');
data_output1=data_output;
load('data_output_figure_4.mat');
data_output4=data_output;

clear data_output;



AlleleReplacement_names={data_output4.strain};

NaturalIsolates_names={data_output1.strain};

for iStrain=1:length(AlleleReplacement_names)
    
    Strain=AlleleReplacement_names{iStrain};
    
    QueryStrain=Strain(6:regexp(Strain,' ')-1);
    
    if strcmp(QueryStrain,'BC')
        
        QueryStrain='BC187';
        
    end
    
    if strcmp(QueryStrain,'YJM')
        
        QueryStrain='YJM978';
        
    end
    
    if strcmp(QueryStrain,'Y9')
        
        QueryStrain='Y9-WashU';
        
    end
    
    idx=find(strcmp(NaturalIsolates_names,QueryStrain));
    
    if idx
        x1(iStrain,1)=mean(data_output4(iStrain).values);
        x1(iStrain,2)=mean(data_output1(idx).values);
    end
end
%%
hfig=figure;
hold all;

x=x1(:,1);
y=x1(:,2);
plot(x,y,'.');


coeffs = polyfit(x, y, 1);
% Get fitted values
fittedX = linspace(min(x), max(x), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 3);


xlabel('Allele replacement set point')
ylabel('Natural isolate set point')
Set_fig_RE(hfig,9,9,20)


[R,P]=corrcoef(x1(:,1),x1(:,2));


Correlation_Coefficient=R(1,2);
P_Value=P(1,2);
