function plot_compute_meiotic_segregants()

% PLOT_COMPUTE_MEIOTIC_SEGREGATNS plots meoitic segregant data
%Inspired by compute_plates_hists_backup in the same path

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
counter_4=1;

% Check that all the fields are the same for all the plates to be combined

plates=fieldnames(plates_bfp);

for a=1:length(plates)
    
    strains=fieldnames(plates_bfp.(plates{a}));
    
%Compute the xi for each Well{iRow,jCol}
%     
%     dat_bfp_yfp=log10(plates_bfp.(plates{a}).(Well{1,12}).FITC_H);
%     [f_norm1,x_norm1] = ksdensity(dat_bfp_yfp);
%     
%     dat_mCh_yfp=log10(plates_mCh.(plates{a}).(Well{1,12}).FITC_H);
%     [f_norm2,x_norm2] = ksdensity(dat_mCh_yfp);
%             
%     dat_other_yfp=log10(plates_other.(plates{a}).(Well{1,12}).FITC_H);
%     [f_norm3,x_norm3] = ksdensity(dat_other_yfp);
    
    %Bins for the calculation of the ks density
    x_norm1=[3.05551957005425,3.06621588075973,3.07691219146521,3.08760850217069,3.09830481287618,3.10900112358166,3.11969743428714,3.13039374499263,3.14109005569811,3.15178636640359,3.16248267710908,3.17317898781456,3.18387529852004,3.19457160922552,3.20526791993101,3.21596423063649,3.22666054134197,3.23735685204746,3.24805316275294,3.25874947345842,3.26944578416390,3.28014209486939,3.29083840557487,3.30153471628035,3.31223102698584,3.32292733769132,3.33362364839680,3.34431995910228,3.35501626980777,3.36571258051325,3.37640889121873,3.38710520192422,3.39780151262970,3.40849782333518,3.41919413404067,3.42989044474615,3.44058675545163,3.45128306615711,3.46197937686260,3.47267568756808,3.48337199827356,3.49406830897905,3.50476461968453,3.51546093039001,3.52615724109549,3.53685355180098,3.54754986250646,3.55824617321194,3.56894248391743,3.57963879462291,3.59033510532839,3.60103141603387,3.61172772673936,3.62242403744484,3.63312034815032,3.64381665885581,3.65451296956129,3.66520928026677,3.67590559097225,3.68660190167774,3.69729821238322,3.70799452308870,3.71869083379419,3.72938714449967,3.74008345520515,3.75077976591064,3.76147607661612,3.77217238732160,3.78286869802708,3.79356500873257,3.80426131943805,3.81495763014353,3.82565394084902,3.83635025155450,3.84704656225998,3.85774287296546,3.86843918367095,3.87913549437643,3.88983180508191,3.90052811578740,3.91122442649288,3.92192073719836,3.93261704790384,3.94331335860933,3.95400966931481,3.96470598002029,3.97540229072578,3.98609860143126,3.99679491213674,4.00749122284223,4.01818753354771,4.02888384425319,4.03958015495867,4.05027646566416,4.06097277636964,4.07166908707512,4.08236539778061,4.09306170848609,4.10375801919157,4.11445432989705];
    x_norm2=[1.39627821431182,1.40505571240787,1.41383321050391,1.42261070859996,1.43138820669600,1.44016570479204,1.44894320288809,1.45772070098413,1.46649819908017,1.47527569717622,1.48405319527226,1.49283069336831,1.50160819146435,1.51038568956039,1.51916318765644,1.52794068575248,1.53671818384852,1.54549568194457,1.55427318004061,1.56305067813666,1.57182817623270,1.58060567432874,1.58938317242479,1.59816067052083,1.60693816861688,1.61571566671292,1.62449316480896,1.63327066290501,1.64204816100105,1.65082565909709,1.65960315719314,1.66838065528918,1.67715815338523,1.68593565148127,1.69471314957731,1.70349064767336,1.71226814576940,1.72104564386544,1.72982314196149,1.73860064005753,1.74737813815358,1.75615563624962,1.76493313434566,1.77371063244171,1.78248813053775,1.79126562863379,1.80004312672984,1.80882062482588,1.81759812292193,1.82637562101797,1.83515311911401,1.84393061721006,1.85270811530610,1.86148561340214,1.87026311149819,1.87904060959423,1.88781810769028,1.89659560578632,1.90537310388236,1.91415060197841,1.92292810007445,1.93170559817049,1.94048309626654,1.94926059436258,1.95803809245863,1.96681559055467,1.97559308865071,1.98437058674676,1.99314808484280,2.00192558293884,2.01070308103489,2.01948057913093,2.02825807722698,2.03703557532302,2.04581307341906,2.05459057151511,2.06336806961115,2.07214556770720,2.08092306580324,2.08970056389928,2.09847806199533,2.10725556009137,2.11603305818741,2.12481055628346,2.13358805437950,2.14236555247555,2.15114305057159,2.15992054866763,2.16869804676368,2.17747554485972,2.18625304295576,2.19503054105181,2.20380803914785,2.21258553724390,2.22136303533994,2.23014053343598,2.23891803153203,2.24769552962807,2.25647302772411,2.26525052582016];
    x_norm3=[0.257417926945458,0.309504009264171,0.361590091582884,0.413676173901597,0.465762256220310,0.517848338539022,0.569934420857735,0.622020503176448,0.674106585495161,0.726192667813874,0.778278750132587,0.830364832451300,0.882450914770013,0.934536997088725,0.986623079407438,1.03870916172615,1.09079524404486,1.14288132636358,1.19496740868229,1.24705349100100,1.29913957331972,1.35122565563843,1.40331173795714,1.45539782027585,1.50748390259457,1.55956998491328,1.61165606723199,1.66374214955071,1.71582823186942,1.76791431418813,1.82000039650684,1.87208647882556,1.92417256114427,1.97625864346298,2.02834472578170,2.08043080810041,2.13251689041912,2.18460297273783,2.23668905505655,2.28877513737526,2.34086121969397,2.39294730201269,2.44503338433140,2.49711946665011,2.54920554896883,2.60129163128754,2.65337771360625,2.70546379592496,2.75754987824368,2.80963596056239,2.86172204288110,2.91380812519982,2.96589420751853,3.01798028983724,3.07006637215595,3.12215245447467,3.17423853679338,3.22632461911209,3.27841070143081,3.33049678374952,3.38258286606823,3.43466894838694,3.48675503070566,3.53884111302437,3.59092719534308,3.64301327766180,3.69509935998051,3.74718544229922,3.79927152461793,3.85135760693665,3.90344368925536,3.95552977157407,4.00761585389279,4.05970193621150,4.11178801853021,4.16387410084892,4.21596018316764,4.26804626548635,4.32013234780506,4.37221843012378,4.42430451244249,4.47639059476120,4.52847667707992,4.58056275939863,4.63264884171734,4.68473492403605,4.73682100635477,4.78890708867348,4.84099317099219,4.89307925331091,4.94516533562962,4.99725141794833,5.04933750026704,5.10142358258576,5.15350966490447,5.20559574722318,5.25768182954190,5.30976791186061,5.36185399417932,5.41394007649803];
    
    
    for iRow=1:8
        for jCol=1:12
            
            if(sum(strcmp(Well(iRow,jCol),strains))==1)
                
                for ichannel=1:length(channels)
                    
                    %plot(plates_hists.(plates{a}).(Well{i,jCol}).(channels{ichannel}).xi,plates_hists.(plates{a}).(Well{i,jCol}).(channels{ichannel}).f)
                    
                    %for each channel compute the corresponding stats
                    
                    % Retrieve data for each channel
                    
                    switch channels{ichannel}
                        
                        case 'bfp_yfp'
                            
                            dat_bfp_yfp=log10(plates_bfp.(plates{a}).(Well{iRow,jCol}).FITC_H);
                            
                            if (length(dat_bfp_yfp)>n_min_events)
                                %BFP

                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).mean=nanmean(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).median=nanmedian(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).perc_ind=sum(dat_bfp_yfp>threshold)./length(dat_bfp_yfp);
                                plates_hists.(plates{a}).(Well{iRow,jCol}).(channels{ichannel}).counts=length(dat_bfp_yfp);
                                
%                                 BC187(counter2)=nanmean(dat_bfp_yfp);
%                                 counter2=counter2+1;
                                
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
                                
%                                 YJM978(counter1)=nanmean(dat_mCh_yfp);
%                                 counter1=counter1+1;
                                
                                
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
                                
                                %Meiotic segregants can be counted on their
                                %own as long as they have good quality data
                                %the parental strains do not matter that
                                %much
                           
                                meioitic_segregants(counter,1)=nanmean(dat_other_yfp);
                                meioitic_segregants(counter,2)=nanmean(dat_bfp_yfp);
                                meioitic_segregants(counter,3)=nanmean(dat_mCh_yfp);
                                
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
                
                %If all strains have quality data
              %if (length(dat_bfp_yfp)>n_min_events & length(dat_mCh_yfp)>n_min_events & length(dat_other_yfp)>n_min_events)
                 
                  meioitic_segregants_4(counter_4,1)=nanmean(dat_other_yfp);
                  meioitic_segregants_4(counter_4,2)=nanmean(dat_bfp_yfp);
                  meioitic_segregants_4(counter_4,3)=nanmean(dat_mCh_yfp);
                  counter_4=counter_4+1;
                  
              %end
                
            end
            
            counts=counts+1;
        end
    end
    
end

%% 20150402 remove cases where mCherry was too high

hfig=figure('Position',[440   562   318   236])

idx_to_remove_2=meioitic_segregants(:,3)>1.85;
meioitic_segregants(idx_to_remove_2,:)=[];

[N,binCenters] = hist(meioitic_segregants(:,1));
hBar = bar(binCenters,N./sum(N),'hist');
set(hBar,'FaceColor',[1,1,1]*0.5,'LineWidth',1) %Fill bars in gray
%xlim([1.5 4.5])
%title(' Distribution of meiotic segregants')

Set_fig_RE(hfig,9,9,20)

YJM978_like=sum(N(1:6))+N(7)./2;
BC187_like=N(7)./2+sum(N(8:end));

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
export_fig(filename,'-pdf',  '-transparent', '-nocrop')

end