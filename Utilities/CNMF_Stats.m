function stats = CNMF_Stats(W,H,X,rmv_flag)

if nargin <4; rmv_flag=0; end

if rmv_flag %Remove empty w or w with barely anything
   [W,H] = RemoveEmptyMotifs(W,H);
end

stats.n_motifs = size(H,1);

Xhat = tensor_convolve(W,H);
Residuals = X-Xhat;

%explained variance
stats.pev = CalculateExplainedVariance(X,Residuals);

%per frame
stats.pev_frame = CalculateExplainedVarianceFrameWise(X,Residuals);

%correlation
stats.rho = corr(X(:),Xhat(:));

%Correlation per frame
stats.rho_frame = zeros(1,size(Xhat,2));
for i = 1:size(Xhat,2)    
    stats.rho_frame(i) = corr(Xhat(:,i),X(:,i));
end    

%loadings of each w
K = size(H,1); 
stats.loadings = zeros(1,K);        
for i = 1:K
    temp = tensor_convolve(W(:,i,:),H(i,:)); 
    stats.loadings(i) = CalculateExplainedVariance(X,X-temp);
end
stats.loadings = stats.loadings/sum(stats.loadings);

%final cost
stats.cost= sqrt(mean((X(:)-Xhat(:)).^2));

%loss
stats.loss = norm(X-Xhat,'fro')/norm(X,'fro');     

%approximate motif frequency

    
end
    
    
    
    
    
    
    
    
    
    
