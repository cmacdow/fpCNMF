function [lambda, fh] = FitLambda(X,W,opts,genfigs)
%Camden MacDowell - timeless. Follows Mackevivious et al., 2018, eLife
if nargin <4; genfigs =1; end

%Sweep Lambda
if opts.verbose; fprintf('\n Fitting Lambda'); end
cost = NaN(1,numel(opts.lambda));
reg = NaN(1,numel(opts.lambda));

if isempty(W)
    for i = 1:numel(opts.lambda)
        if opts.verbose; fprintf('\n\t ........ Fitting Lambda %d of %d ...........', i,numel(opts.lambda)); end
        [W,H] = fpCNMF(X,'L',opts.L,'K',opts.K,'non_penalized_iter',opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
            'speed','fast','verbose',0,'lambda',opts.lambda(i),'ortho_H',opts.ortho_H);
        [cost(i),reg(i)] = CNMF_CostAndReg(X,W,H);
    end
else %just fit Hs
    for i = 1:numel(opts.lambda)
        if opts.verbose; fprintf('\n\t ........ Fitting Lambda %d of %d ...........', i,numel(opts.lambda)); end
        [W,H] = fpCNMF(X,'non_penalized_iter',opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
            'speed','fast','verbose',0,'lambda',opts.lambda(i),'W_update',0,'W',W,'ortho_H',opts.ortho_H);
        [cost(i),reg(i)] = CNMF_CostAndReg(X,W,H);
    end
end
    

%normalize
cost_norm = (cost-min(cost))/(max(cost)-min(cost));
reg_norm = (reg-min(reg))/(max(reg)-min(reg));

%interpolate and find intersection (where crosses zero). could also do this by fitting a polynomial
xq = linspace(opts.lambda(1),opts.lambda(end),1000); 
cost_norm_int = interp1(opts.lambda,cost_norm,xq,'linear');
reg_norm_int = interp1(opts.lambda,reg_norm,xq,'linear');

%approximate intersection
idx = find(reg_norm_int-cost_norm_int<0,1,'first');
lambda = xq(idx);

%optionally make figures;
if genfigs
   fp = fig_params;
   figure; hold on; 
   plot(opts.lambda,cost_norm,'ro','linewidth',1.25); 
   plot(opts.lambda,reg_norm,'bo','linewidth',1.25); 
   set(gca,'xscale','log');
   plot(xq,cost_norm_int,'r','linewidth',1.25);plot(xq,reg_norm_int,'b','linewidth',1.25); plot(lambda,reg_norm_int(idx),'kx','linewidth',2)
   xlabel('Lambda'); ylabel('AU');
   title('Automated Lambda Selection','Fontweight',fp.font_weight,'Fontsize',fp.font_size)
   fh = gcf;
else
   fh = [];
end

end
