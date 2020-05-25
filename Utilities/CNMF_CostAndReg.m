function [cost, reg] = CNMF_CostAndReg(X,W,H)
%from seqNMF toolbox    
[N,K,L] = size(W);
[~,T] = size(H);
kernel = 0.01*ones(1,(2*L)-1);
WTX = zeros(K,T);
for cur_L = 1:L
   WTX = WTX + W(:,:,cur_L)' * [X(:,cur_L:T) zeros(N,cur_L-1)];
end 

reg = (conv2(WTX, kernel, 'same')*H') .* ~eye(K);
reg = norm(reg(:),1);

Xest = tensor_convolve(W,H);
cost = sqrt(sum(X(:)-Xest(:)) .^2);

end %function
