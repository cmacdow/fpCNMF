function [W,H] = CNMF_pMU(X,W,H,terms,update_W)
%update w (N x K x L tensor) and h (K x T matrix) using a penalized multiplicative update. Performs a
%single iteration.

if nargin <5; update_W = 1; end
[N, K, L] = size(W);
[~, T] = size(H);

kernel = 0.01*ones(1,(2*L)-1);  %smoothing kernel

%avoid zero-locking
H = max(eps(),H); 
W = max(eps(),W);

% Compute Estimate
Xest = tensor_convolve(W,H);


% Compute CNMF Update 
WTX = zeros(K,T);
WTXest = zeros(K,T);
for cur_L = 1:L
   WTX = WTX + W(:,:,cur_L)' * circshift(X,1-cur_L,2);
   WTXest = WTXest + W(:,:,cur_L)' * circshift(Xest,1-cur_L,2);
end 

% Update H
if terms.lambda || terms.ortho_H %Optional penalty: X_ortho penalty and H_ortho penalty. 
    penalty = terms.lambda .* ~eye(K) * conv2(WTX, kernel, 'same') + terms.ortho_H * ~eye(K) * conv2(H, kernel, 'same')+ terms.sparse_H;    
else
    penalty = terms.sparse_H; 
end
H = H .* WTX ./ (WTXest + penalty + eps);

%normalize W and H 
[W,H] = NormWH(W,H,N,K,L);

if update_W ==1 %option to skip w_update
    Xest = tensor_convolve(W,H); % Compute Estimate
    if terms.lambda ||terms.ortho_W %precompute constants
        Wflat = sum(W,3);
        XS = conv2(X,kernel,'same');
    end
    for cur_L = 1:L
        H_shift = circshift(H,cur_L-1,2);
        XHT = X * H_shift';
        XestHT = Xest * H_shift';    

        if terms.lambda || terms.ortho_W %Optional penalty: X_ortho penalty and H_ortho penalty. 
            penalty = terms.lambda .* XS * H_shift' * ~eye(K) + terms.ortho_W * Wflat * ~eye(K) + terms.sparse_W;
        else
            penalty = terms.sparse_W;
        end

        W(:,:,cur_L) = W(:,:,cur_L) .* XHT ./ (XestHT + penalty + eps);
    end
end

end

