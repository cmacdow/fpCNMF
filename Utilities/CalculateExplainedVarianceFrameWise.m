function [ExpVar_frame] = CalculateExplainedVarianceFrameWise(X,Residuals)
    %Calculate the variance
    ExpVar_frame = zeros(1,size(Residuals,2));
    for i = 1:size(Residuals,2)    
        ExpVar_frame(i) = 100*(1 - nanvar(Residuals(:,i))/nanvar(X(:,i)));
    end
end