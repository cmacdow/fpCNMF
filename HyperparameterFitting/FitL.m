function [L, fh] = FitL(X,opts,genfigs)
%Camden MacDowell - timeless
if nargin <3; genfigs =1; end

if opts.verbose; fprintf('\n Fitting L values using %d non-penalized iterations',opts.paramsweep_non_penalized_iter); end

%confirm that only a single K exists (otherwise fit L and K together)
assert(numel(opts.K)==1, 'Multiple opts.K values. Either provide a single K value or use FitKandL for combined fit');

%sweep K
stats = cell(1,numel(opts.L));
for i = 1:numel(opts.L)
    if opts.verbose; fprintf('\n\t ........ Fitting L %d of %d ...........', i,numel(opts.L)); end
    [~,~,stats{i}] = fpCNMF(X,'L',opts.L(i),'K',opts.K,'non_penalized_iter',opts.paramsweep_non_penalized_iter,'penalized_iter',0,...
        'speed','fast','verbose',0,'lambda',0);        
end
pev = arrayfun(@(x) stats{x}.pev, 1:numel(stats),'UniformOutput',1);

%Select a point well above the elbow
idx = ceil(1.1*Elbow_pt(pev,[],1));

L = opts.L(idx);

if genfigs
   fp = fig_params;
   figure; hold on; 
   plot(opts.L,pev,'-o','color','r','linewidth',1.25); 
   plot(L,pev(idx),'kx','linewidth',2)
   xlabel('K'); ylabel('PEV');
   title('Automated K Selection','Fontweight',fp.font_weight,'Fontsize',fp.font_size)
   fh = gcf;
else
   fh = [];   
end


end %function end


