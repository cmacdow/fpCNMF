function [lambda, fh] = FitLambda(X,W,opts,genfigs)
%Camden MacDowell - timeless. Follows Mackevivious et al., 2018, eLife
if nargin <4; genfigs =1; end

%Sweep Lambda
if opts.verbose; fprintf('\n Fitting Lambda'); end
cost = NaN(1,numel(opts.lambda_range));
reg = NaN(1,numel(opts.lambda_range));

if isempty(W)
    for i = 1:numel(opts.lambda_range)
        if opts.verbose; fprintf('\n\t ........ Fitting Lambda %d of %d ...........', i,numel(opts.lambda_range)); end
        [W,H] = fpCNMF(X,'L',opts.L,'K',opts.K,'non_penalized_iter',opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
            'speed','fast','verbose',0,'lambda',opts.lambda_range(i));
        [cost(i),reg(i)] = CNMF_CostAndReg(X,W,H);
    end
else %just fit Hs
    for i = 1:numel(opts.lambda_range)
        if opts.verbose; fprintf('\n\t ........ Fitting Lambda %d of %d ...........', i,numel(opts.lambda_range)); end
        [W,H] = fpCNMF(X,'non_penalized_iter',opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
            'speed','fast','verbose',0,'lambda',opts.lambda_range(i),'W_update',0,'W',W);
        [cost(i),reg(i)] = CNMF_CostAndReg(X,W,H);
    end
end
    

%normalize
cost_norm = (cost-min(cost))/(max(cost)-min(cost));
reg_norm = (reg-min(reg))/(max(reg)-min(reg));

%interpolate and find intersection (where crosses zero). could also do this by fitting a polynomial
xq = linspace(opts.lambda_range(1),opts.lambda_range(end),1000); 
cost_norm_int = interp1(opts.lambda_range,cost_norm,xq,'linear');
reg_norm_int = interp1(opts.lambda_range,reg_norm,xq,'linear');

%approximate intersection
idx = find(reg_norm_int-cost_norm_int<0,1,'first');
lambda = xq(idx);

%optionally make figures;
if genfigs
   figure; hold on; 
   plot(opts.lambda_range,cost_norm,'ro'); 
   plot(opts.lambda_range,reg_norm,'bo'); 
   set(gca,'xscale','log');
   plot(xq,cost_norm_int,'r');plot(xq,reg_norm_int,'b'); plot(lambda,reg_norm_int(idx),'kx')
   xlabel('Lambda'); ylabel('AU');
   fh = gcf;
else
   fh = [];
end

end
