function [K, fh] = FitK(X,opts,genfigs)
%Camden MacDowell - timeless
if nargin <3; genfigs =1; end

if opts.verbose; fprintf('\n Fitting K values using %d non-penalized iterations',opts.paramsweep_non_penalized_iter); end

%confirm that only a single L exists (otherwise fit L and K together)
assert(numel(opts.L)==1, 'Multiple opts.L values. Either provide a single L value or use FitKandL for combined fit');

%sweep K
stats = cell(1,numel(opts.K));
for i = 1:numel(opts.K)
    if opts.verbose; fprintf('\n\t ........ Fitting K %d of %d ...........', i,numel(opts.K)); end
    [~,~,stats{i}] = fpCNMF(X,'L',opts.L,'K',opts.K(i),'non_penalized_iter',opts.paramsweep_non_penalized_iter,'penalized_iter',0,...
        'speed','fast','verbose',0,'lambda',0);        
end
pev = arrayfun(@(x) stats{x}.pev, 1:numel(stats),'UniformOutput',1);

%Select a point well above the elbow
idx = ceil(1.1*Elbow_pt(pev,[],1));

K = opts.K(idx);

if genfigs
   fp = fig_params;
   figure; hold on; 
   plot(opts.K,pev,'-o','color','r','linewidth',1.25); 
   plot(K,pev(idx),'kx','linewidth',2)
   xlabel('K'); ylabel('PEV');
   title('Automated K Selection','Fontweight',fp.font_weight,'Fontsize',fp.font_size)
   fh = gcf;
else
   fh = [];   
end


end %function end


