function [W, H, stats_train, stats_test] = FitandValidateMotifs(data_train,data_test,opts,genfigs)
%Camden MacDowell - timeless

%for reproducibility
rng('default');

%Determine the number of non-penalized iterations
if isempty(opts.non_penalized_iter) 
    if opts.verbose; fprintf('\nFitting non penalized iterations'); end
    opts.non_penalized_iter = FitNonPenalizedIterations(data_train,[],opts,genfigs);
end

%Optionally Fit Lambda Penalty Term
if numel(opts.lambda)>1
    if opts.verbose; fprintf('\nFitting Lambda'); end
    lambda = FitLambda(data_train,[],opts,genfigs);
else
    lambda = opts.lambda;
end

%Fit Motifs To Training Data And Collect Statistics
W_temp = cell(1,opts.repeat_fits);
H_temp = cell(1,opts.repeat_fits);
stats_train_temp =cell(1,opts.repeat_fits);
for cur_fit = 1:opts.repeat_fits %fit multiple times due to random initialization
    if opts.verbose; fprintf('\nFitting Training Data Round %d of %d',cur_fit,opts.repeat_fits); end
    [W_temp{cur_fit},H_temp{cur_fit},stats_train_temp{cur_fit}] = fpCNMF(data_train,'L',opts.L,'K',opts.K,'non_penalized_iter',...
        opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
        'speed','fast','verbose',opts.verbose,'lambda',lambda,...
        'ortho_H',opts.ortho_H,'w_update_iter',opts.w_update_iter,...
        'sparse_H',opts.sparse_H);  
    %Remove Empty Motifs 
    [W_temp{cur_fit},H_temp{cur_fit}] = RemoveEmptyMotifs(W_temp{cur_fit},H_temp{cur_fit});
end
%choose best fit
idx = InternallyValidateWs(data_train,W_temp,H_temp,opts.fit_criterion,1);
W = W_temp{idx}; H = H_temp{idx}; stats_train = stats_train_temp{idx}; 
stats_train.lambda = lambda; 


if genfigs %visualize example fit
   VisualizeData(tensor_convolve(W,H),W,H);
   title('Example Fit To Training Data','Fontweight','normal','Fontsize',16); drawnow
end   

%Cross-validate to testing data by only fitting temporal weightings
if opts.verbose; fprintf('\nFitting Testing Data'); end
[Wxval,Hxval,stats_test] = fpCNMF(data_test,'non_penalized_iter',...
    opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
    'speed','fast','verbose',0,'lambda',lambda,'w_update_iter',0,...
    'ortho_H',opts.ortho_H,'W',W,'sparse_H',opts.sparse_H);


if genfigs %visualize example fit
    VisualizeData(tensor_convolve(Wxval,Hxval),Wxval,Hxval);
    title('Example Fit To Testing Data','Fontweight','normal','Fontsize',16); drawnow
end  

end
















