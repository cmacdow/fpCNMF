function [kval,fh] = FitPhenoK(X,varargin)
%Camden MacDowell - timeless. Randomly resamples (with replacement)
%n x n similarity matrix X to determine the best K value for phenographs K
%nearest neighbor search. 

opts.k_range = [2,4,6,8,10,15,20,30,50];
opts.louvain_restarts =5;
opts.verbose = 1;
opts.number_resamples = 25;
opts.genfigs = 1; 
opts = ParseOptionalInputs(opts,varargin);

%for reproducibility
rng('default');

%catch to only use k values at smaller than 1/4th the the size of the data
opts.k_range(opts.k_range>floor(size(X,1))/2)=[]; 

ovr_q = NaN(opts.number_resamples,numel(opts.k_range));
num_clust = NaN(opts.number_resamples,numel(opts.k_range));
for i = 1:numel(opts.k_range) 
    if opts.verbose; fprintf('\n Working on k value %d of %d', i, numel(opts.k_range)); end
    for j = 1:opts.number_resamples
        y = datasample(1:1:size(X,1),size(X,1),'Replace',true);
        [idx, ~, ovr_q(j,i)] = PhenoCluster(X(y,y),'k',opts.k_range(i),'louvain_restarts',opts.louvain_restarts,'Verbose',0);   
        if size(idx,2)>1
            idx = idx(:,end); %phenograph likes to over cluster so use the least clustered level if multiple
        end
        num_clust(j,i) = numel(unique(idx));
    end
end

%Approximate Clusters by taking the elbow in the number of discovered clusters. 
avg_n = nanmean(num_clust);
idx = Elbow_pt(avg_n,[],1);
kval = opts.k_range(idx);

%optionally make figures;
if opts.genfigs
   fp=fig_params;
   figure; hold on; 
   plot(opts.k_range,avg_n,'o-','color','b','linewidth',1.25);
   plot(opts.k_range(idx),avg_n(idx),'X','color','k','linewidth',1.25,'markersize',15);
   ylabel('Number of Clusters');
   xlabel('Number of Neighbors');
   fp.FormatAxes(gca);
   fh = gcf;
else
   fh = [];
end

