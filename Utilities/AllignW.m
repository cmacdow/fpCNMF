function W_alligned = AllignW(W,core_comm_idx,lags,clust_id,lag_mat)

[N, K, L] = size(W);

%pad X
X_pad = NaN(N,K,L*3);
for i =1:K
    X_pad(:,i,:) = cat(2, zeros(N,L), squeeze(W(:,i,:)), zeros(N,L));     %Zero pad
end
W_alligned = X_pad;
for i = 1:numel(unique(clust_id))
    template_idx = core_comm_idx{i}(randperm(numel(core_comm_idx{i}),1)); %get a template pattern from that cluster; 
    cluster_idx = find(clust_id==i);    
    
    %allign all to the template use the lags
    for j = 1:numel(cluster_idx)
        shift_val=lag_mat(cluster_idx(j),template_idx);        
        W_alligned(:,cluster_idx(j),:) = circshift(squeeze(X_pad(:,cluster_idx(j),:)),lags(shift_val),2);
    end
end
    
end