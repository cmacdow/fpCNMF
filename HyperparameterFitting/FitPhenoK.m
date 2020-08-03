function [kval,fh] = FitPhenoK(X,varargin)
%Camden MacDowell - timeless. 

opts.k_range = [2,5,10,15,20,30,50];
opts.louvain_restarts = 3;
opts.verbose = 1;
opts.genfigs = 1; 
opts = ParseOptionalInputs(opts,varargin);

%for reproducibility
rng('default');

%catch to only use k values at smaller than 1/4th the the size of the data
opts.k_range(opts.k_range>floor(size(X,1))/2)=[]; 

num_clust = NaN(1,numel(opts.k_range));
ovr_q = NaN(1,numel(opts.k_range));
tightness = cell(1,numel(opts.k_range));
for i = 1:numel(opts.k_range) 
    if opts.verbose; fprintf('\n Working on k value %d of %d', i, numel(opts.k_range)); end
    [idx,~,ovr_q(i)] = PhenoCluster(X,'k',opts.k_range(i),'louvain_restarts',opts.louvain_restarts,'verbose',1);   
    if size(idx,2)>1; idx = idx(:,end); end %phenograph can return a hierachical cluster levels. this is not validated. don't use.     
    num_clust(i) = numel(unique(idx));
    
    clust = unique(idx);
    for j = 1:numel(clust)
        temp = X(idx==clust(j),idx==clust(j));    
        temp(1:1+size(temp,1):end) = NaN; %remove the autocorrelation
        tightness{i}(j) = nanmean(temp(:));
    end    
end

%Approximate Clusters by taking the elbow in the number of discovered clusters. 
% idx = Elbow_pt(num_clust,[],1);

%Approcimate the choice of K by taking 1 step before the drop off in within cluster correlation
avg_tightness = cellfun(@nanmean,tightness,'UniformOutput',1);
idx = Elbow_pt(avg_tightness,[],1)-1;

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
   plot(opts.k_range,avg_tightness,'o-','color','b','linewidth',1.25);
   ylabel('Cluster Autocorrelation');
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

