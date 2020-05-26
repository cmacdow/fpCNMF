function opts = SyntheticDataOptions(varargin)
%default options for synthetic data fitting example

%General CNMF options
opts.paramsweep_non_penalized_iter = 10;
opts.K = [1,2,3,4,6,10];
opts.L = [5, 10, 15, 20, 25];
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

%Clustering Parameters
opts.clust_smooth_kernel = [];
opts.clust_knn = 2:1:15;
opts.clust_removepad = 1;
opts.clust_num_resamples = 1;
opts.clust_louvain_restarts = 1; 
opts.clust_community_fraction = 1; 

ParseOptionalInputs(opts,varargin);

end