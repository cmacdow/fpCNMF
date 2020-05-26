function FittingResults(stats_train,stats_test,opts)

fp = fig_params; 

%Stats structure may be different sizes depending on fpCNMF speed used. 
%So just take the final statistics from each fit
stats_train = cellfun(@(x) x(1,end),stats_train,'UniformOutput',0);
stats_test = cellfun(@(x) x(1,end),stats_test,'UniformOutput',0);

%concatenate
stats_train = [stats_train{:}];
stats_test = [stats_test{:}];

%compare PEV 
figure; hold on; 
%train
y = [stats_train(:).pev]*100;
x = ones(1,numel(y))*0.9+rand(1,numel(y))*0.2;
plot(x,y,'linestyle','none','marker','.','markersize',15,'color',[0.5 0.5 0.5]);
%test
y = [stats_test(:).pev]*100;
x = 1+ones(1,numel(y))*0.9+rand(1,numel(y))*0.2;
plot(x,y,'linestyle','none','marker','.','markersize',15,'color',[0.5 0 0]);
xlim([0 3]); ylim([50 100]);
ylabel('Percent Explained Variance'); 
set(gca,'Xtick',[1 2],'XtickLabel',{'Train','Test'});
fp.FormatAxes(gca);
fp.SetTitle(gca,'Comparing Explained Variance');
fp.FigureSizing(gcf,[2 2 5 8],[])

%plot the number of discovered motifs
if numel(opts.K)==1
    figure; hold on; 
    h = histogram([stats_train(:).n_motifs],'BinWidth',1);
    maxcount = max(h.Values);
    plot([opts.K,opts.K],[0,maxcount],'linewidth',2,'color','k');
    xlabel('Number of Motifs'); 
    ylabel('Number of Blocks');
    fp.FormatAxes(gca);
    fp.SetTitle(gca,{'Evaluating Number of';'Discovered Motifs'});
    fp.FigureSizing(gcf,[2 2 5 7],[]);
end


end %function end


