function [kval,fh] = FitPhenoK(X,varargin)
%Camden MacDowell - timeless. 

opts.k_range = [2,4,6,8,10,15,20,30,50];
opts.louvain_restarts =5;
opts.verbose = 1;
opts.genfigs = 1; 
opts = ParseOptionalInputs(opts,varargin);

%for reproducibility
rng('default');

%catch to only use k values at smaller than 1/4th the the size of the data
opts.k_range(opts.k_range>floor(size(X,1))/2)=[]; 

num_clust = NaN(1,numel(opts.k_range));
ovr_q = NaN(1,numel(opts.k_range));
for i = 1:numel(opts.k_range) 
    if opts.verbose; fprintf('\n Working on k value %d of %d', i, numel(opts.k_range)); end
    [idx,~,ovr_q(i)] = PhenoCluster(X,'k',opts.k_range(i),'louvain_restarts',opts.louvain_restarts,'Verbose',0);   
    if size(idx,2)>1; idx = idx(:,end); end %phenograph can return a hierachical cluster levels. don't use.     
    num_clust(i) = numel(unique(idx));
end

%Approximate Clusters by taking the elbow in the number of discovered clusters. 
idx = Elbow_pt(num_clust,[],1);

%if multiple knn values produce the same number of clusters as the elbow, then use the one with the highest modularity
same_values = find(num_clust == num_clust(idx));
if numel(same_values)>1
   [~, best_q] = max(ovr_q(same_values));
   idx = same_values(best_q(1));
end

kval = opts.k_range(idx);

%optionally make figures;
if opts.genfigs
   fp=fig_params;
   figure; hold on; 
   yyaxis left
   plot(opts.k_range,num_clust,'o-','color','k','linewidth',1.25);
   plot(opts.k_range(idx),num_clust(idx),'X','color','k','linewidth',1.25,'markersize',15);
   ylabel('Number of Clusters');
   yyaxis right
   plot(opts.k_range,ovr_q,'o-','color','b','linewidth',1.25);
   ylabel('Louvain Q Modularity');
   xlabel('Number of Neighbors');
   title('Autofitting PhenoCluster KNN','FontSize',fp.font_size,'Fontweight',fp.font_weight);
   ax = gca;
   ax.YAxis(1).Color = 'k';
   ax.YAxis(2).Color = 'b';
   fp.FormatAxes(gca);
   fh = gcf;
else
   fh = [];
end

