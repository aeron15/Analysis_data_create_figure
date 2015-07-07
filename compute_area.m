function [perc_area] = compute_area(y,x,off_dist_y,off_dist_x);

% The COMPUTE_METRICS function is designed to compute fraction of induced
% cells, based off of a predefined 'OFF' peak--which is designated as a
% glucose only lane

% Inputs: for query = y (ks-density.f), x (ks-density.xi)

% modified by KL 20150615 to work specifically for single-gradient data 

% set up a common x 

if isnan(off_dist_x)|isinf(off_dist_x)
    perc_area = nan;
    
else
    max_x = max([max(x),max(off_dist_x)]);
    min_x = min([min(x),min(off_dist_x)]);
    xx = linspace(min_x,max_x,150);
    % 
    try
        yy = interp1(x,y,xx,'linear');yy(yy<0)=0;
        yy(isnan(yy))=0;
        yy_off = interp1(off_dist_x,off_dist_y,xx);
        yy_off(yy_off<0)=0;
        yy_off(isnan(yy_off))=0;
        % xx = off_dist_x;

        % normalizes the two (query and reference) distributions to an area that
        % equals 1
        pp = yy/trapz(xx,yy);
        p_off =yy_off/trapz(xx,yy_off);
        diff_p = pp-p_off;

        [val,ind] = max(p_off);
        % xx(ind)
        diff_p(diff_p<0)=0;
        diff_p(1:ind) = 0;

        % plot(xx,p_off);hold on;
        % plot(xx,diff_p,'r';
        perc_area = trapz(xx,diff_p);
        % plot(xx,p_off,'r',xx,pp,'k',xx,diff_p)
    catch
        perc_area = nan;
    end
end