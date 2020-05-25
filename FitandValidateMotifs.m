function [W, stats_train, stats_test] = FitandValidateMotifs(data_train,data_test,opts,genfigs)
%Camden MacDowell - timeless

%Fit K and L (pending)

%Determine the number of non-penalized iterations
opts.non_penalized_iter = DetermineNonPenalizedIterations(data_train,[],opts,genfigs);

%Fit Lambda
if genfigs %generate figures for an example block
    lambda = FitLambda(data_train{block},[],opts,1);
    title('Automated Lambda Selection','Fontweight','normal','Fontsize',16);
else %no figures
    lambda = FitLambda(data_train{block},opts,0);
end

%Fit Motifs To Training Data And Collect Statistics
[W{block},H,stats_train{block}] = fpCNMF(data_train,'L',opts.L,'K',opts.K,'non_penalized_iter',...
    opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
    'speed','fast','verbose',0,'lambda',lambda);

%Remove Empty Motifs 
[W{block},H] = RemoveEmptyMotifs(W,H);

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
    sgtitle('Example Fit To Testing Data, Block','Fontweight','normal','Fontsize',16); drawnow
end  

end
















