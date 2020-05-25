function [W_basis, kval, ovr_q, idx_louvain, idx_knn, tcorr_mat, handles] = ClusterW(W,opts)
%camden macdowell - timeless
%file_list is the full path of each file with motifs to cluster. motifs are
%contained in variable 'w' as a pxl x motif x time tensor

if iscell(W) %combine mutliple CNMF Fits stored in a cell array
   W = cat(2,W{:});
end

%Remove empty motifs 
W = RemoveEmptyMotifs(W);

[N, ~, L] = size(W);

%optional 2D or 3D gaussian smooth. Reccomended for noisy data. 
if ~isempty(opts.clust_smooth_kernel)
    W_smooth = GaussianSmoothTensor(W,opts.clust_smooth_kernel);
else
    W_smooth = W; 
end

%compute maximum temporal xcorrelation between motifs
[tcorr_mat, lag_mat, lags] = TemporalXCorrTensor(W_smooth,L,1);

%If given a range of neighbors fit Phenograph Number of Neighbors
if numel(opts.clust_knn)>1
    [kval, ~] = FitPhenoK(tcorr_mat,'k_range',opts.clust_knn,'louvain_restarts',opts.clust_louvain_restarts,'genfigs',1,...
        'verbose',1,'number_resamples',opts.clust_num_resamples);
else
    kval = opts.clust_knn;
end

%Cluster
[idx_louvain, idx_knn, ovr_q] = PhenoCluster(tcorr_mat,'k',kval,'louvain_restarts',opts.clust_louvain_restarts,'Verbose',1);   
if size(idx_louvain,2) > 1; idx_louvain = idx_louvain(:,end); end %just take the initital clustering. Phenograph tends to really overcluster at the higher levels. 

%visualize cross correlation matrix
Plot_OrderedSimilarityMatrix(tcorr_mat,idx_louvain);

%Get core community to average for basis motifs
[core_comm_idx, ~] = CoreCommunity(idx_louvain,idx_knn,opts.clust_community_fraction); 

%Allign motifs in each cluster to one of the core community members 
W_alligned = AllignW(W,core_comm_idx,lags,idx_louvain,lag_mat);

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


















        