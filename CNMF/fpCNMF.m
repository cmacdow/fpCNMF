function [W, H, stats] = fpCNMF(X,varargin)

%% Initialize
%parse optional inputs
opts.W = []; 
opts.K = 4;
opts.L =15;
opts.non_penalized_iter = 25;
opts.w_update_iter = 1;
opts.speed = 'fast';
opts.penalized_iter = 25;
opts.verbose = 0; 

%specific terms for pMU
opts.lambda = 0;
opts.ortho_H = 0;
opts.ortho_W = 0;
opts.sparse_H = 0;
opts.sparse_W = 0;

opts = ParseOptionalInputs(opts,varargin);

[N, T] = size(X);

%initialize W, if exists, then adjust K and L accordingly
if isempty(opts.W)
    W = max(X(:))*rand(N, opts.K, opts.L);
else
    W = opts.W;
    [~, opts.K, opts.L] = size(W);
end

%initialize X. Padding to prevent edge effects
X = [zeros(N,opts.L), X, zeros(N,opts.L)];

%initialize H
H = [zeros(opts.K,opts.L),max(X(:))*rand(opts.K,T),zeros(opts.K,opts.L)];

%normalize
[W,H] = NormWH(W,H);

%generate vector of fit types (could be expanded to include different algs
fit_type = [ones(1,opts.non_penalized_iter),zeros(1,opts.penalized_iter)];
update_w_idx = cat(2,[1:opts.w_update_iter:opts.non_penalized_iter],[1:opts.w_update_iter:opts.penalized_iter]);
update_w = zeros(1,numel(fit_type));
update_w(update_w_idx)=1;

%% CNMF

switch opts.speed 
    case 'fast-gpu' %barebones. just fit. no reconstructions
        if opts.verbose; fprintf('\n\tUNLEASH THE GLORIOUS POWER OF THE GPU! Only final statistics computed\n');  end
        X = gpuArray(X); W = gpuArray(W); H = gpuArray(H);
        for iter = 1:numel(fit_type)                 
            if fit_type(iter) == 1 
                [W,H] = CNMF_ALS(X,W,H,update_w(iter));
            else
                [W,H] = CNMF_pMU(X,W,H,opts,update_w(iter));
            end
        end
        X = gather(X); W = gather(W); H = gather(H);
        stats =CNMF_Stats(W,H,X,1);
    
    case 'fast' %barebones. just fit. no reconstructions
        if opts.verbose; fprintf('\n\tZoom zoom! Speediness activated! Only final statistics computed\n');  end
        for iter = 1:numel(fit_type)        
            if fit_type(iter) == 1 
                [W,H] = CNMF_ALS(X,W,H,update_w(iter));
            else
                [W,H] = CNMF_pMU(X,W,H,opts,update_w(iter));
            end
        end
        stats =CNMF_Stats(W,H,X,1);
        
        
    case 'normal' %compute basic statistics and evaluate loss per iteration.
        if opts.verbose; fprintf('\n\tNormal Speed! Computing statistics per iteration\n');  end           
        stats = CNMF_Stats(W,H,X,1);        
        for iter = 1:numel(fit_type)                 
            if fit_type(iter) == 1 
                [W,H] = CNMF_ALS(X,W,H,update_w(iter));
            else
                [W,H] = CNMF_pMU(X,W,H,opts,update_w(iter));
            end
            stats(1,iter+1) = CNMF_Stats(W,H,X,1);
        end
                
            
    case 'diagnostic'
        if opts.verbose; fprintf('\n\tDiagnostic Speed, generating figures - hope you figure out what is going on.....');  end       
        stats = CNMF_Stats(W,H,X,1);  
        fh=figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
        for iter = 1:numel(fit_type)  
            tic
            if fit_type(iter) == 1 
                [W,H] = CNMF_ALS(X,W,H,update_w(iter));
            else
                [W,H] = CNMF_pMU(X,W,H,opts,update_w(iter));
            end
            stats(iter+1) = CNMF_Stats(W,H,X,1);
            %visualize           
            VisualizeData(tensor_convolve(W,H),W,H,fh);
        end                   
        
        
    otherwise
        error('Unknown speed. Try again.')
end

%remove pad
H = H(:,opts.L+1:end-opts.L);

end %function

