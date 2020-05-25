function Plot_CompareWs(W1,W2)
%Camden MacDowell - timeless. 
%usually make W1 the original, and W2 the new

%confirm the same number of motifs
assert(size(W1,2)==size(W2,2),'Number of Motifs (K) of compared Ws do not match');

%match the duration (if L was longer than original);
max_L = max(size(W1,3),size(W2,3));

[N,K,L] = size(W1); 
if L<max_L; W1 = cat(3,W1,zeros(N,K,max_L-L)); end

[N,K,L] = size(W2); 
if L<max_L; W2 = cat(3,W2,zeros(N,K,max_L-L)); end

L = max_L;

%find the closest motifs in W2 that match W1
rho = NaN(1,K);
idx = NaN(1,K);
for k = 1:K
    tcorr_mat = TemporalXCorrTensor(cat(2,W1(:,k,:),W2),ceil(L/2),0);
    [rho(k),idx(k)] = max(tcorr_mat(1,2:end)); %ignore autocorrelation
end

%gutcheck warning if multiple motifs in W1 match same W2 (can occur with large noise)
if numel(unique(idx))<K; warning('multiple W1s matched to same W2'); end

W1 = ShiftW(W1); W2 = ShiftW(W2);

%reorder W2 to match W1
W2 = W2(:,idx,:);
rho = rho(idx);

figure; hold on; 
ax1 = subplot(1,2,1); hold on; 
ax2 = subplot(1,2,2); hold on; 

axes(ax1); VisualizeW(W1,1,0); title('Original','fontweight','normal')
axes(ax2); VisualizeW(W2,1,0); title([sprintf('Basis Motifs \nRho: '),sprintf('%.2g  ',rho)],'fontweight','normal')


end %function end