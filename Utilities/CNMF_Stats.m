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
temploadings = zeros(1,K);   
stats.rho_frame_per_motif = zeros(size(H,1),size(Xhat,2)); %per each motif (used to compute frequency)
for i = 1:K
    temp = tensor_convolve(W(:,i,:),H(i,:)); 
    temploadings(i) = CalculateExplainedVariance(X,X-temp);
    for j = 1:size(temp,2)    
        stats.rho_frame_per_motif(i,j) = corr(temp(:,j),X(:,j));            
    end        
end
stats.loadings = {temploadings/sum(temploadings)};

%final cost
stats.cost= sqrt(mean((X(:)-Xhat(:)).^2));

%loss
stats.loss = norm(X-Xhat,'fro')/norm(X,'fro');     

    
end
    
    
    
    
    
    
    
    
    
    
