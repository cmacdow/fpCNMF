%% example analysis pipeline for synthetic data

%Generate Synthetic Data
rng('default');
[data, W_orig] = GenerateSyntheticData('verbose',1,'N',50,'T',500,'sigma_data',0.025,'L',15,'K',4,'blocks',10); sgtitle('Example Data Block');

%Split into equal numbers of test/train blocks
data_train = data(1:2:end);
data_test = data(2:2:end);
data_train = data_train(1:min(numel(data_train),numel(data_test)));
data_test = data_test(1:min(numel(data_train),numel(data_test)));

%get configurable fitting options
opts = SyntheticDataOptions();

%fit to all the recording blocks
stats_train = cell(1,numel(data_train));
stats_test = cell(1,numel(data_train));
W = cell(1,numel(data_train));

%get an (over) approximation of K and/or L using 1 example block
if numel(opts.K)>1 && numel(opts.L)==1 %fit K for a given L
    opts.K = FitK(data_train{1},opts,1);
elseif numel(opts.L)>1 && numel(opts.K)==1 %fit L for a given K
    opts.L = FitL(data_train{1},opts,1);
elseif numel(opts.L)>1 && numel(opts.K)>1 %fit together (usually not necessary)
     [opts.K, opts.L] = FitKandL(data_train{1},opts,1);
else %Use the existing K and L
end

%run example fits in parallel 
parpool(parcluster('local'));
parfor block = 1:numel(data_train)
    fprintf('\n Working on block %d of %d',block,numel(data_train));   
    if block ==1 %generate example figures
        [W{block}, ~, stats_train{block}, stats_test{block}] = FitandValidateMotifs(data_train{block},data_test{block},opts,1);    
    else
        [W{block}, ~, stats_train{block}, stats_test{block}] = FitandValidateMotifs(data_train{block},data_test{block},opts,0);    
    end
end

% Visualize General Motifs Statistics
FittingResults(stats_train,stats_test,opts);

% Cluster Motifs
W_basis = ClusterW(W,opts);

% Compare the basis motifs to the original motifs
Plot_CompareWs(W_orig,W_basis);

%Refit to all the data
data_all = cat(2,data{:}); %entire "recording"

%Recompute statistics that can change as a function of numel(data);
opts.non_penalized_iter = FitNonPenalizedIterations(data_all,W_basis,opts,1);

lambda = FitLambda(data_all,W_basis,opts,1); title('Automated Lambda Selection, Entire Recording','Fontweight','normal','Fontsize',16);

%Refit temporal weightings
[~,H,stats] = fpCNMF(data_all,'non_penalized_iter',...
    opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
    'speed','normal','verbose',0,'lambda',lambda,'w_update_iter',0,'W',W_basis);
    
%Visualize
VisualizeData(data_all,W_basis,H); sgtitle('Example Fit To Entire Recording','Fontweight','normal','Fontsize',16); 
CNMF_Stats(W_basis,H,data_all);























