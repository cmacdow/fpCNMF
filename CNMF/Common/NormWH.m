function [W,H] = NormWH(W,H)

norms = sqrt(sum(H.^2, 2))';
H = diag(1 ./ (norms+eps)) * H;
W = permute(reshape(reshape(permute(W,[1 3 2]),[],size(W,2))*diag(norms),size(W,1),size(W,3),[]),[1 3 2]); %faster than for loop for each L
    
end %function

