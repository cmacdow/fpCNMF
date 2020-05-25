function fh = VisualizeData(data,W,H,fh)

if nargin < 4; fh=figure('units','normalized','position',[0.1 0.1 0.8 0.8]); end

set(0, 'currentfigure', fh);
hold on; 

[N,K,~] = size(W);
[~,T] = size(H);

col = distinguishable_colors(K);

%Create subplots
ax1 = subplot(1,3,1); hold on; title('W: Spatiotemporal Motifs','FontWeight','normal')
ax2 = subplot(1,3,2); hold on; title('H: Temporal Weightings','FontWeight','normal')
ax3 = subplot(1,3,3); hold on; title('Reconstruction','FontWeight','normal')

%w
axes(ax1)
VisualizeW(W,1,0);
set(ax1,'units','normalized','position',[0.05 0.1 0.15 0.65])

%h
axes(ax2)
for i = 1:K
    plot(H(i,:),'linewidth',2,'color',col(i,:));
end
set(ax2,'xtick',[],'ytick',[]); 
set(ax2,'units','normalized','position',[0.25 0.80 0.65 0.1])

%data
axes(ax3)
imagesc(data); colormap(flipud(gray));
ylim([0.5,N+0.5]); xlim([0,T]); axis off; box on
set(ax3,'units','normalized','position',[0.25 0.1 0.65 0.65])

%link axes
linkaxes([ax2,ax3],'x');


fh = gcf; 
end


