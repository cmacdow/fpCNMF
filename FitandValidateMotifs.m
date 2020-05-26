function [W, stats_train, stats_test] = FitandValidateMotifs(data_train,data_test,opts,genfigs)
%Camden MacDowell - timeless

%Optionally Fit K and L. Seperately or together
if numel(opts.K)>1 && numel(opts.L)==1 %fit K for a given L
    opts.K = FitK(data_train,opts,genfigs);
elseif numel(opts.L)>1 && numel(opts.K)==1 %fit L for a given K
    opts.L = FitL(data_train,opts,genfigs);
elseif numel(opts.L)>1 && numel(opts.K)>1 %fit together (usually not necessary)
     [opts.K, opts.L] = FitKandL(data_train,opts,genfigs);
else %Use the existing K and L
end

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
















