function fh = VisualizeH(H)
%H is a K x T matrix

[K,~] = size(H);
col = distinguishable_colors(K);
figure; hold on; 
label = cell(1,K);
for i = 1:K
    plot(H(i,:),'linewidth',2,'color',col(i,:));
    label{i} = sprintf('Motif %d',i);
end
legend(label{:});
fh = gcf;

end