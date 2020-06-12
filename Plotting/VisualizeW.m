function fh = VisualizeW(W,concat_flag,new_fig_flag)
%W is a P,K,L tensor
if nargin<2; concat_flag=0; end
if nargin<3; new_fig_flag=1; end
     
[N,K,L] = size(W);
col = distinguishable_colors(K);

%generate new figure (vs plot on existing)
if new_fig_flag; figure('Units','normalize','position',[0.1 0.1 0.8 0.8]); hold on; end

if concat_flag
   gap = 2;
   W_concat = zeros(N,(L+gap)*K); %add small gap between Ws
   idx = cell(1,K);
   for i = 1:K
       temp = cat(2,squeeze(W(:,i,:)),zeros(N,gap));       
       temp = temp/max(temp(:));
       W_concat(:,((i-1)*(L+gap)+1) : (i*(L+gap))) = cat(2,squeeze(W(:,i,:)),zeros(N,gap));       
       idx{i} = [((i-1)*(L+gap))+0.5,(i*(L+gap))-1.5]; %get indices of motif
   end
   W_concat(:,end-gap)=[]; %remove last pad
   imagesc(W_concat); colormap(flipud(gray));    
   %outline colors
   for i = 1:K
       plot(idx{i},[0.5 0.5],'linewidth',4,'color',col(i,:))
       plot([idx{i}(2),idx{i}(2)],[0 N+1],'linewidth',2,'color',col(i,:))
   end
   xlim([0, size(W_concat,2)]);
   ylim([0.5, size(W_concat,1)+0.5]);
   axis off;
   fh = gcf;
else %make into a seperate figure with each motif a subpanel
    [nR,nC] = numSubplot(K,2);    
    for i = 1:K
       subplot(nR,nC,i); 
       imagesc(squeeze(W(:,i,:))); 
       colormap(flipud(gray)); 
       title(sprintf('Motif %d',i),'FontSize',14,'FontWeight','normal');
       axis equal; axis off
    end
    fh = gcf;
end

end
    
   