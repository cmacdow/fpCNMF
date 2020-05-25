function opts = SyntheticDataOptions(varargin)
%default options for synthetic data fitting example

%DiscoveringMotifs_SytheticData()
opts.clust_smooth_kernel = [];
opts.clust_knn = 2:1:15;
opts.clust_removepad = 1;
opts.clust_num_resamples = 10;
opts.clust_louvain_restarts = 1; 
opts.clust_community_fraction = 1; 

%parse optional inputs
opts.K = 6;
opts.L = 20;
opts.max_non_penalized_iter =10; 
opts.w_update_iter = 1;
opts.speed = 'fast';
opts.penalized_iter = 100;
opts.verbose = 0; 

%specific terms for pMU
opts.lambda = sort(logspace(0,-5,15), 'ascend');
opts.ortho_H = 0;
opts.ortho_W = 0;
opts.sparse_H = 0;
opts.sparse_W = 0;

ParseOptionalInputs(opts,varargin);

end