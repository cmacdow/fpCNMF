function [ExpVar_all] = CalculateExplainedVariance(X,Residuals)

%Across all frames
ExpVar_all = 1 - nanvar(Residuals(:))./nanvar(X(:));

end