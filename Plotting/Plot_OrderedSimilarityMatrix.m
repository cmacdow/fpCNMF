function [fh, x_ordered] = Plot_OrderedSimilarityMatrix(x,clust_id,varargin)
figure; hold on; 

%Camden MacDowell - timeless
fp = fig_params;
%parse optional inputs
opts.order = sort(unique(clust_id),'ascend'); %order just by the cluster_id
opts.colormap = 'magma';
opts.range = [0.1 0.9];
opts.offset = 0;
opts.ytextoffset = -1;
opts.xtextoffset = -2;
opts.label_thickness = 2;
opts.fig_position = [680   558   560   420];
opts.label_col = distinguishable_colors(numel(unique(clust_id)));
opts = ParseOptionalInputs(opts,varargin);

%gut check
assert(numel(opts.order)==numel(unique(clust_id)),'number of elements in order does not match number of clusters')

%sort
location = [];
for i = 1:numel(opts.order)
    location = cat(1,location,find(clust_id==opts.order(i)));
end
x_ordered = x(location,location);
clust_id_ordered = clust_id(location);

%visualize
imagesc(x_ordered,opts.range); hold on
set(gca,'Ydir','reverse','xlim',[0.5,size(x_ordered,2)+0.5],'ylim',[0.5,size(x_ordered,1)+0.5])
colormap(opts.colormap)
set(gcf,'Position',opts.fig_position)

tickval = [];
for i = 1:numel(opts.order)
    %find middle coord
    tickval =(find(clust_id_ordered == i));
    tickval = linspace(min(tickval)-0.5,max(tickval)+0.5,numel(tickval));
    %Add colored lines 
    plot(ones(1,numel(tickval))*-1*opts.offset,tickval,'LineWidth',opts.label_thickness,'Color',opts.label_col(i,:))
    %Add X label
    text(opts.offset+opts.ytextoffset,mean(tickval),['\color[rgb]{' num2str(opts.label_col(i,:)),sprintf('}%d',i)],'HorizontalAlignment','Center','FontSize',16)   
    %Add Y label
    text(mean(tickval),numel(clust_id)-opts.xtextoffset+opts.offset,['\color[rgb]{' num2str(opts.label_col(i,:)),sprintf('}%d',i)],'rotation',0,'HorizontalAlignment','Center','FontSize',16)
    %Add colored lines 
    plot(tickval,ones(1,numel(tickval))+numel(clust_id)+opts.offset,'LineWidth',opts.label_thickness,'Color',opts.label_col(i,:))
end

cB = colorbar;
ylabel(cB, 'Similarity','FontSize',16,'FontWeight','normal')
set(cB,'YTick',linspace(opts.range(1),opts.range(2),3))
set(gca,'box','off','Clipping','off','XColor','none','YColor','none')
set(gca,'XTick',tickval,'YTick',tickval,'TickDir','out','XTickLabelRotation',90);
xlabel('Cluster IDs');
ylabel('Cluster IDs');
title('Similarity Between Discovered Motifs','FontWeight',fp.font_weight,'FontSize',fp.font_size);
axis square
fp.FormatAxes(gca);

fh = gcf;


end %function

