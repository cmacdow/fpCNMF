function [W,H] = NormWH(W,H,N,K,L)

norms = sqrt(sum(H.^2, 2))';
H = diag(1 ./ (norms+eps)) * H;
W = permute(reshape(reshape(permute(W,[1 3 2]),[],K)*diag(norms),N,L,[]),[1 3 2]); %faster than for loop for each L

end %function

