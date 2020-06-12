function [K, L, fh] = FitKandL(X,opts,genfigs)
%Camden MacDowell - timeless. Follows Mackevivious et al., 2018, eLife
if nargin <4; genfigs =1; end

%Sweep K and L together
if opts.verbose; fprintf('\n Fitting L and K'); end

sweep = combvec(opts.L,opts.K);

stats = cell(1,size(sweep,2));
for i = 1:size(sweep,2)
    if opts.verbose; fprintf('\n\t ........ Fitting L/K %d of %d ...........', i,size(sweep,2)); end
    [~,~,stats{i}] = fpCNMF(X,'L',sweep(1,i),'K',sweep(2,i),'non_penalized_iter',opts.max_non_penalized_iter,'penalized_iter',0,...
        'speed','fast','verbose',0,'lambda',0,'W_update',opts.w_update_iter);        
end
pev = arrayfun(@(x) stats{x}.pev, 1:numel(stats),'UniformOutput',1);

% get peak pev 
[~,idx] = max(pev);
L = ceil(sweep(1,idx)*1.1);
K = ceil(sweep(2,idx)*1.1);

%optionally make figures;
if genfigs
   fp = fig_params;
   figure;
   [xq,yq] = meshgrid(linspace(min(sweep(1,:)),max(sweep(1,:)),100),linspace(min(sweep(2,:)),max(sweep(2,:)),100));
   vq = griddata(sweep(1,:),sweep(2,:),pev,xq,yq);
   mesh(xq,yq,vq)
   hold on
   plot3(sweep(1,:),sweep(2,:),pev,'o')
   plot3([L,L],[K,K],[0 1],'color','k','linewidth',2)   
   xlabel('L'); ylabel('K'); zlabel('PEV');
   title('Automated K and L Selection','Fontweight',fp.font_weight,'Fontsize',fp.font_size)
   fh = gcf;
else
   fh = [];
end

end
