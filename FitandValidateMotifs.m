function [W, stats_train, stats_test] = FitandValidateMotifs(data_train,data_test,opts,genfigs)
%Camden MacDowell - timeless

%Determine the number of non-penalized iterations
opts.non_penalized_iter = FitNonPenalizedIterations(data_train,[],opts,genfigs);

%Optionally Fit Lambda Penalty Term
if numel(opts.lambda)>1
    lambda = FitLambda(data_train,[],opts,genfigs);
end

%Fit Motifs To Training Data And Collect Statistics
[W,H,stats_train] = fpCNMF(data_train,'L',opts.L,'K',opts.K,'non_penalized_iter',...
    opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
    'speed','fast','verbose',0,'lambda',lambda);

%Remove Empty Motifs 
[W,H] = RemoveEmptyMotifs(W,H);

if genfigs %visualize example fit
    VisualizeData(data_train,W,H); 
    sgtitle('Example Fit To Training Data','Fontweight','normal','Fontsize',16); drawnow
end   

%Cross-validate to testing data by only fitting temporal weightings
[Wxval,Hxval,stats_test] = fpCNMF(data_test,'non_penalized_iter',...
    opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
    'speed','fast','verbose',0,'lambda',lambda,'w_update_iter',0,'W',W);

if genfigs %visualize example fit
    VisualizeData(data_test,Wxval,Hxval);
    sgtitle('Example Fit To Testing Data','Fontweight','normal','Fontsize',16); drawnow
end  

end
















