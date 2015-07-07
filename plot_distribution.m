function plot_distribution(x)

hfig=figure('Position',[440   562   318   236])
[N,binCenters] = hist(x);
hBar = bar(binCenters,N./sum(N),'hist');
set(hBar,'FaceColor',[1,1,1]*0.5,'LineWidth',1) %Fill bars in gray
Set_fig_RE(hfig,9,9,20)

YJM978_like=sum(N(1:6))+N(7)./2;
BC187_like=N(7)./2+sum(N(8:end));

end