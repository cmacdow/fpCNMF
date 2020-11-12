function [core_comm_idx,core_comm_size,idx,crit_val] = AutoFitCommunityNumber(cluster_idx,idx_knn,opts,W_smooth,lag_mat,lags)
    criteria = NaN(numel(unique(cluster_idx)),numel(opts.clust_community_fraction));
    for cur_lvl = 1:numel(opts.clust_community_fraction) 
        [core_comm_idx, ~] = CoreCommunity(cluster_idx,idx_knn,opts.clust_community_fraction(cur_lvl),1); 

        %Allign motifs in each cluster to one of the core community members 
        W_alligned = AllignW(W_smooth,core_comm_idx,lags,cluster_idx,lag_mat);

        %compute basis motifs
        basis_temp = NaN(size(W_smooth,1),numel(core_comm_idx),size(W_alligned,3));
        for i = 1:numel(core_comm_idx)   
            basis_temp(:,i,:) = nanmean(W_alligned(:,core_comm_idx{i},:),2);
        end
        good_idx = nanvar(squeeze(sum(basis_temp,1)),[],1)>eps;
        basis_temp = basis_temp(:,:,good_idx); 
                
        switch opts.autofit_method
            case 'autocorrelation' %compute autocorrelation to each motif in cluster
                for cur_basis = 1:size(basis_temp,2)
                    motif_idx = find(cluster_idx==cur_basis);
                    flat_motifs = NaN(numel(basis_temp(:,1,:)),numel(motif_idx));
                    for cur_motif = 1:numel(motif_idx)
                        temp = W_alligned(:,motif_idx(cur_motif),good_idx);
                        flat_motifs(:,cur_motif) = temp(:);
                    end
                    temp = basis_temp(:,cur_basis,:);
                    criteria(cur_basis,cur_lvl) = nanmean(corr(temp(:),flat_motifs));
                end
                critfun = @(x) max(x,[],2);
            case 'stvar' %spatiotemporal variance
                for cur_basis = 1:size(basis_temp,2)
                    motif_idx = find(cluster_idx==cur_basis);
                    flat_motifs = NaN(numel(basis_temp(:,1,:)),numel(motif_idx));
                    for cur_motif = 1:numel(motif_idx)
                        temp = W_alligned(:,motif_idx(cur_motif),good_idx);
                        flat_motifs(:,cur_motif) = temp(:);
                    end
                    temp = basis_temp(:,cur_basis,:);
                    criteria(cur_basis,cur_lvl) = abs(nanmean(nanvar(temp(:))-nanvar(flat_motifs,[],1)));
                end    
                critfun = @(x) min(x,[],2);
        end
    end %the the community level loop 
    [crit_val,idx] = critfun(criteria);
    
    %now compute each core community separately
    core_comm_idx = cell(1,numel(idx));
    core_comm_size = NaN(1,numel(idx));
    for cur = 1:numel(idx)
        %get the dyanmicsness of the average motifs
        [temp, ~] = CoreCommunity(cluster_idx,idx_knn,opts.clust_community_fraction(idx(cur)),1); 
        core_comm_idx{cur} = temp{cur};
        core_comm_size(cur) = opts.clust_community_fraction(idx(cur));
    end
    
end