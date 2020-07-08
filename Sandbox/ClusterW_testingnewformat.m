function [W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs] = ClusterW_testingnewformat(W,opts,nanpxs)
%camden macdowell - timeless
%file_list is the full path of each file with motifs to cluster. motifs are
%contained in variable 'w' as a pxl x motif x time tensor

%optional 2D or 3D gaussian smooth. Reccomended. Required if using different Nanpixels across animals.  
W_smooth = W; %The current use of both W_smooth and W is combersome, but it allows different smoothing kernels to be used for clustering vs averaging. 
if iscell(nanpxs) && iscell(W) && ~isempty(opts.clust_smooth_kernel)  %if nanpxs is a cell, then you need to smooth data individual    
    fprintf('\t\n Smoothing Data using individual nan masks\n')
    for i = 1:numel(nanpxs)
        W_smooth{i} = GaussianSmoothTensor(W{i},opts.clust_smooth_kernel,opts.originaldimensions,nanpxs{i});
        temp_smooth = zeros(opts.originaldimensions(1)*opts.originaldimensions(2),size(W_smooth{i},2),size(W_smooth{i},3));
        temp = zeros(opts.originaldimensions(1)*opts.originaldimensions(2),size(W_smooth{i},2),size(W_smooth{i},3));
        for k = 1:size(W_smooth{i},2) %recondition individual cells of W so same size
            temp_smooth(:,k,:) = reshape(conditionDffMat(squeeze(W_smooth{i}(:,k,:))',nanpxs{i},[],[opts.originaldimensions,size(W_smooth{i}(:,k,:),3)]),[size(temp_smooth,1),1,size(temp_smooth,3)]);
            temp(:,k,:) = reshape(conditionDffMat(squeeze(W{i}(:,k,:))',nanpxs{i},[],[opts.originaldimensions,size(W_smooth{i}(:,k,:),3)]),[size(temp_smooth,1),1,size(temp_smooth,3)]);
        end
        W_smooth{i} = temp_smooth;
        W{i} = temp; 
    end
elseif ~iscell(nanpxs) && iscell(W) && ~isempty(opts.clust_smooth_kernel) %same Nan for all 
    fprintf('\t\n Smoothing Data using same nan mask (cell)\n')
    for i = 1:numel(nanpxs)
        W_smooth{i} = GaussianSmoothTensor(W{i},opts.clust_smooth_kernel,opts.originaldimensions,nanpxs);
    end
elseif ~iscell(nanpxs) && ~iscell(W) && ~isempty(opts.clust_smooth_kernel) %no cell loop
    fprintf('\t\n Smoothing Data using same nan mask (cell)\n')
     W_smooth = GaussianSmoothTensor(W,opts.clust_smooth_kernel,opts.originaldimensions,nanpxs);
else %do nothing
end

if iscell(W) %catch for about if is in cell form. Combine mutliple CNMF Fits stored in a cell array
   W = cat(2,W{:});
   W_smooth = cat(2,W_smooth{:});
end

[N, ~, L] = size(W);  

%recompute the nanpxs across all data using the smoothed data
W_smooth(isnan(W_smooth))=0;
W(isnan(W))=0;
nanpxs = find(nanvar(reshape(W_smooth,[size(W,1),size(W,2)*size(W,3)]),[],2)<=eps);
W_smooth(nanpxs,:,:) = [];
W(nanpxs,:,:) = [];

%compute maximum temporal xcorrelation between motifs
[tcorr_mat, lag_mat, lags] = TemporalXCorrTensor_BigData(W_smooth,ceil(L/2),1);

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


















        