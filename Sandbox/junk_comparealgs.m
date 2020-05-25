function [W, H] = fpCNMF(X,varargin)

%parse optional inputs
opts.K = 4;
opts.L =15;
opts.non_penalized_iter = 25;
opts.w_update_iter = 4;
opts.speedy = 1;
opts.penalized_iter = 25;
opts.verbose = 1; 

opts = ParseOptionalInputs(opts,varargin);

[N, T] = size(X);

%initialize and pad h (randomized and normalized) and w (using ANLS); 
X = [zeros(N,opts.L), X, zeros(N,opts.L)];
H = [zeros(opts.K,opts.L),max(X(:))*rand(opts.K,T)./(sqrt(T/3)),zeros(opts.K,opts.L)]; 
W = max(X(:))*rand(N, opts.K, opts.L); 
[W,H] = NormWH(W,H,N,opts.K,opts.L);

H_init = H;
W_init = W;

% 
% if opts.speedy
%     fprintf('\n\tZoom zoom! Speediness activated! No plots, Iteration Counting, or Iteration statistics\n');          
% end

W = W_init;
H = H_init;
ppev = [];
aloss = [];
acost = [];
stats = CNMFStats(W,H,X,0);
Xhat = tensor_convolve(W,H);
acost(1) = sqrt(mean((X(:)-Xhat(:)).^2));
aloss(1) = stats.loss;
at = 0;
terms = GetpMUTerms(opts.L);
for iter = 1:opts.penalized_iter 
    tic
    if iter<=4
        [W,H] = CNMF_ALS(X,W,H,1); %mod(iter,opts.w_update_iter)+1        
    else
        [W,H] = CNMF_pMU(X,W,H,terms,1);
    end
    stats = CNMFStats(W,H,X,0);
    Xhat = tensor_convolve(W,H);
    acost(iter+1) = sqrt(mean((X(:)-Xhat(:)).^2));
    aloss(iter+1) = stats.loss;
    at(iter+1) = toc;
    if sum(at)>30
        break
    end
end
% 
% figure; hold on; title('loss');  plot(aloss,'k'); set(gca,'yscale','log');
% figure; hold on; title('cost');  plot(acost,'k');set(gca,'yscale','log');


W = W_init;
H = H_init;
pev = [];
loss = [];
cost = [];
stats = CNMFStats(W,H,X,0);
Xhat = tensor_convolve(W,H);
cost(1) = sqrt(mean((X(:)-Xhat(:)).^2));
pev(1) = stats.pev;
loss(1) = stats.loss;
t=0;
terms = GetpMUTerms(opts.L);
for iter = 1:opts.penalized_iter 
    tic
    [W,H] = CNMF_pMU(X,W,H,terms,1);
    stats = CNMFStats(W,H,X,0);
    Xhat = tensor_convolve(W,H);
    cost(iter+1) = sqrt(mean((X(:)-Xhat(:)).^2));
    loss(iter+1) = stats.loss;
%     VisualizeData(Xhat,W,H); pause(); close;
    t(iter+1)=toc;
    if sum(t)>30
        break
    end
end

figure; hold on; title('loss'); plot(cumsum(t),loss,'r'); plot(cumsum(at),aloss,'k'); set(gca,'yscale','log');legend('MU','ANLS');
figure; hold on; title('cost'); plot(cumsum(t),cost,'r'); plot(cumsum(at),acost,'k');set(gca,'yscale','log');legend('MU','ANLS');


figure; hold on; title('loss'); plot(loss,'r'); plot(aloss,'k'); set(gca,'yscale','log');legend('MU','ANLS');
figure; hold on; title('cost'); plot(cost,'r'); plot(acost,'k');set(gca,'yscale','log');legend('MU','ANLS');


%remove zeropadding
% H = H(:,opts.L+1:end-opts.L);

end %function end
% 
%     
% end
% 
% 
% tic
% L = opts.L;
% K = opts.K;
% %Run CNMF_ALS + CNMF_pMU
% for iter = 1:20 %% add the padding
%     if iter ==1
%         X_pad = [zeros(N,L),X,zeros(N,L)];
%         w = w_init; 
%         h = [zeros(K,L),h_init,zeros(K,L)];
%         Xhat = helper.reconstruct(w,h);       
%         cost = sqrt(mean((X_pad(:)-Xhat(:)).^2));
%         [w,h] = CNMF_ALS(X_pad,w,h);   
%     else
%         [w,h] = CNMF_ALS(X_pad,w,h);   
%     end
%        
%     Xhat = helper.reconstruct(w,h);   
%     cost(iter+1) = sqrt(mean((X_pad(:)-Xhat(:)).^2));
% %     loss(iter+1) = norm(X_pad-Xhat,'fro')/norm(X,'fro');
%     %fit until the loss plateuas    
% end
% h = h(:,L+1:end-L);
% t = toc;
% 
% %switch the CNMF_pMU 
% [w, h, temp] = seqNMF(X, ...    % X is the data matrix
%   'K', opts.K, 'L', opts.L, 'lambda', 0, ...        % Other inputs optional
%   'showPlot', 0, 'maxiter', 10, 'tolerance', 0, 'shift', 0, ... 
%   'lambdaL1W', 0, 'lambdaL1H', 0, ...
%   'lambdaOrthoH', 0, 'lambdaOrthoW', 0,...
%   'W_init',w,'H_init',h);
% temp = temp(2:end);
% 
% 
% %Run CNMF_MU_Penalized
% 
% %switch the CNMF_pMU 
% tic
% [w, h, cost2] = seqNMF(X, ...    % X is the data matrix
%   'K', opts.K, 'L', opts.L, 'lambda', 0, ...        % Other inputs optional
%   'showPlot', 0, 'maxiter', 20, 'tolerance', 0, 'shift', 0, ... 
%   'lambdaL1W', 0, 'lambdaL1H', 0, ...
%   'W_init', w_init,'H_init', h_init,...
%   'lambdaOrthoH', 0, 'lambdaOrthoW', 0);
% cost2 = cost2(1:end);
% t2 = toc;
%   
% 
% figure; cost = [cost, temp'];
% plot(cost); hold on;
% plot(cost2); hold on;
% 
% 
% 
% figure;hold on; plot((0:20)*t1,cost,'k'); plot((0:20)*t2,cost2,'b');set(gca,'yscale','log'); ylim([0 3])
% plot(cost); hold on;
% plot(cost2); hold on;

