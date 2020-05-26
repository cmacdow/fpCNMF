function epsilon = FitDBSCANepsilon(D, minpts, genfigs)
%follows kNNdistplot() from dbscan python package

if nargin <3; genfigs = 1; end

%compute the average distance of every point ot it's k nearest neighbors
distance = sort(nanmean(maxk(D,minpts,2),2),'ascend'); 

%get the elbow
idx = Elbow_pt(distance);
epsilon = distance(idx);

%optionally make figures;
if genfigs
   fp=fig_params;
   figure; hold on; 
   plot(distance,'o-','color','b','linewidth',1.25);
   plot(idx,distance(idx),'X','color','k','linewidth',1.25,'markersize',15);
   ylabel(sprintf('minpoints=%d distance',minpts));
   xlabel('Sorted Points');
   title('Autofitting DBSCAN epsilon','FontSize',fp.font_size,'Fontweight',fp.font_weight);
   fp.FormatAxes(gca);
   fh = gcf;
else
   fh = [];
end


end %function