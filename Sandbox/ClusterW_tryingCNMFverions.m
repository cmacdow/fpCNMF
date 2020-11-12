function [W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs] = ClusterW(W,opts,nanpxs)
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
    W_smooth = GaussianSmoothTensor(W,opts.clust_smooth_kernel,opts.originaldimensions,nanpxs,opts.clust_nobleed);
else
    W_smooth = W; 
end

%Reconstruct full W. This is used for averaging to make the basis motifs. 
W = MaskTensor(W,nanpxs,[opts.originaldimensions(1)*opts.originaldimensions(2),size(W,2),size(W,3)]); 

%compute maximum temporal xcorrelation between motifs
[tcorr_mat, lag_mat, lags] = TemporalXCorrTensor_BigData(W_smooth,ceil(L/2),1);

%Cluster using desired method
switch opts.clust_method
    case 'PhenoCluster' %use phenograph
        %compute maximum temporal xcorrelation between motifs
        [tcorr_mat, lag_mat, lags] = TemporalXCorrTensor_BigData(W_smooth,ceil(L/2),1);
        if numel(opts.clust_knn)>1 %optionally estimate a reasonable k value
            [kval, ~] = FitPhenoK(tcorr_mat,'k_range',opts.clust_knn,'louvain_restarts',opts.clust_louvain_restarts,'genfigs',1,'verbose',0);
        else; kval = opts.clust_knn; 
        end
        [cluster_idx, idx_knn, ovr_q] = PhenoCluster(tcorr_mat,'k',kval,'louvain_restarts',opts.clust_louvain_restarts,'verbose',0);   
        if size(cluster_idx,2) > 1; cluster_idx = cluster_idx(:,end); end %just take the initital clustering level 
        
    case 'DBSCAN' %use DBSCAN
        diss_tcorr_mat = 1-tcorr_mat; %DBSCAN uses a dissimilarity matrix. 
        if isempty(opts.clust_epsilon)
           opts.clust_epsilon = FitDBSCANepsilon(diss_tcorr_mat, opts.clust_minpts, 1);
        end
        cluster_idx = dbscan(diss_tcorr_mat,opts.clust_epsilon,opts.clust_minpts,'Distance','precomputed');
        
    case 'CNMF' %In progress. Uses a second round of fpCNMF to cluster motifs. 
        temp = arrayfun(@(n) cat(2,squeeze(W_smooth(:,n,:)),zeros(size(W_smooth,1),size(W_smooth,3))),1:size(W_smooth,2),'UniformOutput',0);
        temp = cat(2,temp{:});
        
        %Fit Lambda Penalty Term
        lambda = FitLambda(data_train,[],opts,genfigs);


        %Fit Motifs To Training Data And Collect Statistics
        W_temp = cell(1,opts.repeat_fits);
        H_temp = cell(1,opts.repeat_fits);
        stats_train_temp =cell(1,opts.repeat_fits);
        for cur_fit = 1:opts.repeat_fits %fit multiple times due to random initialization
            if opts.verbose; fprintf('\nFitting Training Data Round %d of %d',cur_fit,opts.repeat_fits); end
            [W_temp,H_temp,stats_train_temp] = fpCNMF(temp,'L',13,'K',20,'non_penalized_iter',...
                opts.non_penalized_iter,'penalized_iter',50,...
                'speed','fast','verbose',opts.verbose,'lambda',0.001,...
                'ortho_H',opts.ortho_H,'w_update_iter',opts.w_update_iter,...
                'sparse_H',opts.sparse_H);  
            %Remove Empty Motifs 
            [W_temp{cur_fit},H_temp{cur_fit}] = RemoveEmptyMotifs(W_temp{cur_fit},H_temp{cur_fit});
        end
        %choose best fit
        idx = InternallyValidateWs(data_train,W_temp,H_temp,opts.fit_criterion,1);
        W = W_temp{idx}; H = H_temp{idx}; stats_train = stats_train_temp{idx}; 
        stats_train.lambda = lambda;     
        
        
    otherwise
        error('%s is an unsupported clustering method',opts.clust_method);
end
        
        
%visualize cross correlation matrix
fprintf('\n\tPlotting similarity matrix')
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);

fprintf('\n\tGenerating Basis Motifs')
%Get core community to average for basis motifs
if numel(opts.clust_community_fraction)>1 %find the best core_community_fraction
    fprintf('\n\t Autofitting Community Fractions');
    [core_comm_idx, ~] = AutoFitCommunityFraction(cluster_idx,idx_knn,opts,W_smooth,lag_mat,lags);        
else %just take the set value
    fprintf('\n\t Using Set Community Fraction');
    [core_comm_idx, ~] = CoreCommunity(cluster_idx,idx_knn,opts.clust_community_fraction); 
end

%Allign motifs in each cluster to one of the core community members 
fprintf('\n\tAlligning W')
W_alligned = AllignW(W,core_comm_idx,lags,cluster_idx,lag_mat);

%compute basis motifs
fprintf('\n\tComputing Basis Motifs')
W_basis = NaN(size(W,1),numel(core_comm_idx),size(W_alligned,3));
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


















        