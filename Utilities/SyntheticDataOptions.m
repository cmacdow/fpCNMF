function opts = SyntheticDataOptions(varargin)
%default options for synthetic data fitting example

%general options
opts.originaldimensions = 1; %the pixel x pixel dimensions of each frame. set to 1 if not working with imaging data

%General CNMF options
opts.paramsweep_non_penalized_iter = 10;
opts.K = [1,2,3,4,6,10];
opts.L = 20;
opts.non_penalized_iter = [];
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

%General clustering Parameters
opts.clust_method = 'PhenoCluster';
opts.clust_smooth_kernel = [];
opts.clust_community_fraction = 1; 
opts.clust_removepad = 1;

%PhenoCluster parameters
opts.clust_knn = 2:1:15;
opts.clust_louvain_restarts = 5; 

%DBSCAN parameters
opts.clust_epsilon = 0.3; 
opts.clust_minpts = 4; 

ParseOptionalInputs(opts,varargin);

end