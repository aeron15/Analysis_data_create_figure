function compute_range_of_variation()

%COMPUTE_RANGE_OF_VARIATION computes the variation in the set point of
%induction. THis responds to fold-change.

load data_output_figure_1.mat

%%
BC187_idx=find(strcmp({data_output.strain},'BC187'));
YJM978_idx=find(strcmp({data_output.strain},'YJM978'));

2.^abs(mean(data_output(BC187_idx).values)-mean(data_output(YJM978_idx).values))

BC187_mean_fig1=mean(data_output(BC187_idx).values);
YJM978_mean_fig1=mean(data_output(YJM978_idx).values);


BC187_std_fig1=std(data_output(BC187_idx).values);
YJM978_std_fig1=std(data_output(YJM978_idx).values);


BC187_sem_fig1=std(data_output(BC187_idx).values)./sqrt(length(data_output(BC187_idx).values))
YJM978_sem_fig1=std(data_output(YJM978_idx).values)./sqrt(length(data_output(YJM978_idx).values))


%%
BC187_YJM978_error_difference_Strains=sqrt((BC187_std_fig1).^2+(YJM978_std_fig1).^2);
%Range of the set points on figure 1
YJM421_last_loc_range=2.^abs(mean(data_output(loc(1)).values)-mean(data_output(loc(end)).values));

%Range of the strains on figure 4
YJM421_DBPVG1373_range=2.^abs(mean(data_output(27).values)-mean(data_output(6).values));


%Range of the strains on figure 5

CLIB382_idx=find(strcmp({data_output.strain},'CLIB382'));


YJM978_CLIB382_range=2.^abs(mean(data_output(CLIB382_idx).values)-mean(data_output(YJM978_idx).values));
