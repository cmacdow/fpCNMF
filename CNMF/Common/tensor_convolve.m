function X = tensor_convolve(W,H)
%from helper.reconstruct function in seqNMF toolbox

[N,K,L] = size(W);
[~,T] = size(H);

% zeropad by L
H = [zeros(K,L),H,zeros(K,L)];
T = T+2*L;
X = zeros(N,T);

for tau = 1:L 
    X = X + W(:, :, tau) * circshift(H,[0,tau-1]);
end

% undo zer0padding
X = X(:,(L+1):(end-L)); 
