function [W, H, stats_train, stats_test] = NestedCNMF(path_processed,opts,genfigs)
%Camden MacDowell - timeless

%for reproducibility
rng('default');

%grab the fit motifs for this group
[folder_st,file_st] = fileparts(path_processed);
[fn_motifs,~] = GrabFiles([file_st,'\w*chunk\w*.mat'],0,{folder_st}); 

%load the motifs 
w = {};
for i =1:numel(fn_motifs)
    temp = load(fn_motifs{i},'w');
    w{i} = temp.w;
end
w = cat(2,w{:});

%normalize
for i = 1:size(w,2)
    temp = squeeze(w(:,i,:));
    w(:,i,:) = (temp-min(temp(:)))./(max(temp(:))-min(temp(:)));
end

%get a subset of motifs pad with zeros and refit
n_motifs = 40;
n_runs = 10;
for i = 1:n_runs %subsample 10 random collections of motifs
   idx = randperm(size(w,2),n_motifs);
   temp = arrayfun(@(n) cat(2,squeeze(w(:,n,:)),zeros(size(w,1),15)), idx,'UniformOutput',0);
   temp = cat(2,temp{:});
   [W_temp{i}, H{i}, stats_train{i}, stats_test{i}] = FitIter(temp,opts,0);
   W{i} = cat(2,W_temp{i});
end


end %function 


function [W, H, stats_train, stats_test] = FitIter(data_train,opts,genfigs)
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

end
















