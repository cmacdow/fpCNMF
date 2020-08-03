function [idx, idx_knn, ovr_q]= PhenoCluster(X,varargin)
%adapted from Tim Buschman Phenograph. X is n x n a similarity matrix

%parse optional inputs
opts.k = 5;
opts.louvain_restarts =5;
opts.verbose = 0;
opts = ParseOptionalInputs(opts,varargin);

%compute the knn indices
idx_knn = nan(size(X,1),opts.k);
for i = 1:size(X,1)
   [~,idx_knn(i,:)] = maxk(X(i,:),opts.k); 
end

% Now apply Jacard similarity to create connectivity graph
if opts.verbose; fprintf('Calculating weights of network graph...\n'); end
sel_N = size(idx_knn, 1);
w = sparse(sel_N, sel_N);
for i = 1:sel_N
    %Find those indices with overlap
    intersect_idx = sum(ismember(idx_knn, idx_knn(i, :)), 2);
    j_ind = (intersect_idx > 0);
    w(i, j_ind) = intersect_idx(j_ind)./(2*opts.k - intersect_idx(j_ind));
    if mod(i, round(sel_N/10)) == 0 && opts.verbose, fprintf('\t%d%% done...\n', round(i/sel_N*100)); end
end %node loop
w = (w+w') - eye(size(w,1)).*diag(w); %symmetrize 
if opts.verbose; fprintf('done.\n'); end

% Louvain community detection to find clusters
if opts.verbose; fprintf('Clustering...'); end
[idx, ovr_q] = LouvainCommunity(w, 'NumRandomStarts', opts.louvain_restarts,'Verbose',opts.verbose);
if opts.verbose; fprintf('done.\n'); end

for i = 1:size(idx,2)
    if opts.verbose; fprintf('Level %d: Found %d communities...\n', i, length(unique(idx(:, size(idx, 2) - i + 1)))); end
end

end