function [W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles] = ClusterW(W,opts,nanpxs)
%camden macdowell - timeless
%file_list is the full path of each file with motifs to cluster. motifs are
%contained in variable 'w' as a pxl x motif x time tensor

if nargin <3; nanpxs = []; end %masked pixels for imaging data

if iscell(W) %combine mutliple CNMF Fits stored in a cell array
   W = cat(2,W{:});
end

%Remove empty motifs 
W = RemoveEmptyMotifs(W);

[N, ~, L] = size(W);

%optional 2D or 3D gaussian smooth. Reccomended for noisy data. 
if ~isempty(opts.clust_smooth_kernel)
    W_smooth = GaussianSmoothTensor(W,opts.clust_smooth_kernel,opts.originaldimensions,nanpxs);
else
    W_smooth = W; 
end

%compute maximum temporal xcorrelation between motifs
[tcorr_mat, lag_mat, lags] = TemporalXCorrTensor(W_smooth,L,1);

%Cluster using desired method
switch opts.clust_method
    case 'PhenoCluster' %use phenograph
        if numel(opts.clust_knn)>1 %optionally estimate a reasonable k value
            [kval, ~] = FitPhenoK(tcorr_mat,'k_range',opts.clust_knn,'louvain_restarts',opts.clust_louvain_restarts,'genfigs',1,'verbose',1);
        else; kval = opts.clust_knn; 
        end
        [cluster_idx, idx_knn, ovr_q] = PhenoCluster(tcorr_mat,'k',kval,'louvain_restarts',opts.clust_louvain_restarts,'Verbose',1);   
        if size(cluster_idx,2) > 1; cluster_idx = cluster_idx(:,end); end %just take the initital clustering level 
        
    case 'DBSCAN' %use DBSCAN
        diss_tcorr_mat = 1-tcorr_mat; %DBSCAN uses a dissimilarity matrix. 
        if isempty(opts.clust_epsilon)
           opts.clust_epsilon = FitDBSCANepsilon(diss_tcorr_mat, opts.clust_minpts, 1);
        end
        cluster_idx = dbscan(diss_tcorr_mat,opts.clust_epsilon,opts.clust_minpts,'Distance','precomputed');
    otherwise
        error('%s is an unsupported clustering method',opts.clust_method);
end
        
        
%visualize cross correlation matrix
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);

%Get core community to average for basis motifs
[core_comm_idx, ~] = CoreCommunity(cluster_idx,idx_knn,opts.clust_community_fraction); 

%Allign motifs in each cluster to one of the core community members 
W_alligned = AllignW(W,core_comm_idx,lags,cluster_idx,lag_mat);

%compute basis motifs
W_basis = NaN(N,numel(core_comm_idx),size(W_alligned,3));
for i = 1:numel(core_comm_idx)   
    W_basis(:,i,:) = nanmean(W_alligned(:,core_comm_idx{i},:),2);
end

%optional removal of the padded pixels
if opts.clust_removepad
   W_basis = ShiftW(W_basis); %center shift before optional 
   W_basis = W_basis(:,:,L+1:L*2);
else %just remove totally empty regions
   W_basis = W_basis(:,:,nanvar(squeeze(sum(W_basis,1)),[],1)>eps);
end

handles = get(groot, 'Children');

end


















        