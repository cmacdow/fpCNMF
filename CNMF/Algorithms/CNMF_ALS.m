function [W,H] = CNMF_ALS(X,W,H,update_W)
%update w (N x K x L tensor) and h (K x T matrix) using ALS. Performs a
%single iteration.

if nargin <4; update_W = 1; end

[N, K, L] = size(W);
[~, T] = size(H);

%update H
for t = 1:T 
    last = min(t+L-1,T);
    block_size = last - t + 1;  
    
    b = X(:,t:last);
    b = b(:);
    W_blocked = ConditionW(W,[N*L,K]);
    H(:,t) = nnlsm_blockpivot(W_blocked(1:block_size*N,:), b,0);  
end

%normalize H and W
[W,H] = NormWH(W,H);

%update W
h_block = BlockH(H,L); %make H into a block matrix format

if update_W ==1 %option to skip w_update if set to 0; 
    %Compute the nnls
    W_stacked = nnlsm_blockpivot(h_block', X',0); 
    %stack w
    W = ConditionW(W_stacked',[N,K,L]); 
end

end

