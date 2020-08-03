function [idx, crit] = InternallyValidateWs(X,Ws,Hs,criterion,verbose)
%Camden MacDowell -timeless
%Ws and Hs are cell arrays from multiple fits


N = numel(Ws);
crit = zeros(1,N);
for i = 1:N   %loop through fits
    k = size(Ws{i},2);
    
    switch criterion
        case 'AIC' %Residuals are normally distributed so we can approximate AIC=nlog(?^2Z)+2k,
            Xhat = tensor_convolve(Ws{i},Hs{i});        
            r_var = nanvar(X(:)-Xhat(:));
            crit(i) = numel(X(:))*log(r_var)+(2*k); 
            targetfun = @min;
            %         figure; histogram(X(:)-Xhat(:));
        case 'BIC'
            Xhat = tensor_convolve(Ws{i},Hs{i});          
            r_var = nanvar(X(:)-Xhat(:));
            crit(i) = k*log(numel(X)) - (2*log(r_var));
            targetfun = @min;
        case 'PEV'
            crit(i) = CalculateExplainedVariance(X,X-tensor_convolve(Ws{i},Hs{i}));  
            targetfun = @max;
        otherwise 
            error('unknown crossvalidation criterion'); 
    end %switch

end %fit loop

[~, idx] = targetfun(crit); %get the index with the best value 

%get the best fit (lowest 
if verbose 
    figure; hold on; plot(crit,'linewidth',2,'color',[0.5 0.5 0.5],'marker','o'); 
    plot(idx,crit(idx),'marker','x','linestyle','none','color','r','markersize',20,'linewidth',2);
    ylabel(sprintf('%s (AU)',criterion)); xlabel('fit iteration');     
end

end %function
