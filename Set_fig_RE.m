function []=Set_fig_RE(fig,axis_fs,title_fs,label_fs)

%SET_FIGURE_RE formats figures

% Set figure parameters
set(gcf, 'color', [1 1 1]) ;
N = get(fig,'children');


for i = 1:length(N);
    set(N(i),'fontsize',axis_fs,...
        'fontname','Helvetica');
    
    x = get(N(i),'xlabel');set(x,'fontsize',label_fs,'fontname','Helvetica')
    y = get(N(i),'ylabel');set(y,'fontsize',label_fs,'fontname','Helvetica')
    t = get(N(i),'Title');set(t,'fontsize',title_fs,'fontname','Helvetica')
    
    
end