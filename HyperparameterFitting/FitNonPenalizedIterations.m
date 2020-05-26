function [num_iter, fh]= FitNonPenalizedIterations(X,W,opts,verbose)
%Camden MacDowell - timeless

if nargin <4; verbose = 1; end

if isempty(W)
    [~,~,stats] = fpCNMF(X,'L',opts.L,'K',opts.K,'non_penalized_iter',...
        opts.max_non_penalized_iter,'penalized_iter',0,...
        'speed','normal','verbose',0,'lambda',0);
else
    [~,~,stats] = fpCNMF(X,'non_penalized_iter',...
        opts.max_non_penalized_iter,'penalized_iter',0,...
        'speed','normal','verbose',0,'lambda',0,'w_update_iter',0,'W',W);
end

cost = [stats(:).cost];

%get the elbow
num_iter = Elbow_pt(cost(2:end),[],1)+1; %add 1 to account for initial cost measurement

if verbose
    fp = fig_params;
    figure; hold on; 
    plot(1:numel(cost),cost,'o-','linewidth',1.25,'color','r');
    %mark selected point. 
    plot(num_iter,cost(num_iter),'x','linestyle','none','markersize',10,'color','k','linewidth',2);
    xlabel('Fit Iteration');
    ylabel('Cost');
    title('Selected Number of Non-Penalized Iterations','FontWeight','normal');
    set(gca,'yscale','log');
    fp.FormatAxes(gca);
else
    fh = [];
end



end

