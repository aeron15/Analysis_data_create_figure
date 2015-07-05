function meioitic_segregants=plot_compute_meiotic_segregants()

% PLOT_COMPUTE_MEIOTIC_SEGREGATNS plots meoitic segregant data
% returns counts of BC187 like and YJM978. Computes mean and area metric.

path_data='/Users/RenanEscalante/Dropbox/Phenotypic_diversity/var_facs/20120427_meiotic_segregant_classification/output/';

load([path_data 'plates_bfp.mat']);
load([path_data 'plates_mCh.mat']);
load([path_data 'plates_other.mat']);

%Dependency on a map of 96 well plates
load('map_plate_96');

channels={'bfp_yfp','mCh_yfp','other_yfp'};

%% Pick the threshold of induction for percentage calculation

threshold=2.5;
n_min_events=10;
counts=1;
counter=1;

%% Import off peak distribution GET OFF DISTRIBUTION
offStrain_data=plates_other.('Plate_Plate_1').('C03').FITC_H;
offStrain_data=plates_other.('Plate_Plate_10').('B07').FITC_H;


[offStrain_y1,offStrain_x1]=ksdensity(log10(offStrain_data));
% figure;
% plot(offStrain_x1,offStrain_y1);

%%

% Check that all the fields are the same for all the plates to be combined

plates=fieldnames(plates_bfp);

for a=1:length(plates)
    
    strains=fieldnames(plates_bfp.(plates{a}));
    for iRow=1:8
        for jCol=1:12
            
            if(sum(strcmp(Well(iRow,jCol),strains))==1)
                
                for ichannel=1:length(channels)
                    
                    switch channels{ichannel}
                        
                        case 'bfp_yfp'
                            
                            dat_bfp_yfp=log10(plates_bfp.(plates{a}).(Well{iRow,jCol}).FITC_H);
                            
                            if (length(dat_bfp_yfp)>n_min_events)
                                %BFP
                                
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_bfp_yfp>threshold)./length(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=length(dat_bfp_yfp);
                                
                            else
                                
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_bfp_yfp>threshold)./length(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=0;
                            end
                            
                        case 'mCh_yfp'
                            
                            dat_mCh_yfp=log10(plates_mCh.(plates{a}).(Well{iRow,jCol}).FITC_H);
                            
                            if (length(dat_mCh_yfp)>n_min_events)
                                
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_mCh_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_mCh_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_mCh_yfp>threshold)./length(dat_mCh_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=length(dat_mCh_yfp);
                                
                            else
                                
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).f=0;
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).xi=0;
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_mCh_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_mCh_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_mCh_yfp>threshold)./length(dat_mCh_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=0;
                                
                            end
                            
                        case 'other_yfp'
                            
                            dat_other_yfp=log10(plates_other.(plates{a}).(Well{iRow,jCol}).FITC_H);
                            
                            if (length(dat_other_yfp)>n_min_events)
                                
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_other_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_other_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_other_yfp>threshold)./length(dat_other_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=length(dat_other_yfp);
                                
                                meioitic_segregants(counter,1)=nanmean(dat_other_yfp);
                                meioitic_segregants(counter,2)=nanmean(dat_bfp_yfp);
                                meioitic_segregants(counter,3)=nanmean(dat_mCh_yfp);
                                
                                
                                other_strain=dat_other_yfp;
                                other_strain(isnan(other_strain))=[];
                                
                                if mean(other_strain)<1.5
                                   display('ajlksjalk') 
                                end
                                    
                                   
                                
                                [y,x]=ksdensity(other_strain);
                                [perc_area] = compute_area(y,x,offStrain_y1,offStrain_x1);
                                meioticSegregants_area(counter,1)=perc_area;
                                
                                counter=counter+1;
                                
                            else
                                
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).f=0;
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).xi=0;
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_other_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_other_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_other_yfp>threshold)./length(dat_other_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=0;
                                
                            end
                            
                    end
                    
                    
                end
            end
            counts=counts+1;
        end
    end
    
end

%% Remove cases where mCherry was too high

idx_to_remove_2=meioitic_segregants(:,3)>1.85;

meioitic_segregants(idx_to_remove_2,:)=[];
meioticSegregants_area(idx_to_remove_2,:)=[];

meioticSegregants_mean=meioitic_segregants(:,1);

%%
plot_distribution(meioticSegregants_mean)
%%
plot_distribution(meioticSegregants_area)

%% The area metric is much more biased
figure;
%plot(meioticSegregants_mean,meioticSegregants_area,'.')

scatterhist(meioticSegregants_mean,meioticSegregants_area)
axis square;
filename=['Distribution of meiotic segregants_mean_vs_area_more_events.pdf'];
export_fig_specific_path(filename,'-pdf',  '-transparent', '-nocrop')

%% Compute chi-square statistic

E1=sum(N)./2;
O1=YJM978_like;
%O1=450;

E2=sum(N)./2;
O2=BC187_like;
%O2=906-O1;

chi_square=(O1-E1)^2./E1+(O2-E2)^2./E2;
p=1-chi2cdf(chi_square,1)

%%
filename=['Distribution of meiotic segregants.pdf'];
export_fig_specific_path(filename,'-pdf',  '-transparent', '-nocrop')

end