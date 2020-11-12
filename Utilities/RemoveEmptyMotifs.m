function [W,H] = RemoveEmptyMotifs(W,H)

if nargin <2; H =[]; end

indempty = sum(sum(W>eps,1),3)==0; % W is literally empty
indempty = indempty | (max(sum(W,3),[],1).^2> .5*sum(sum(W,3).^2,1)); % or one pixel has >50% of the power

temp = logical(size(indempty));
for i = 1:size(W,2)
   temp(i) = (max(squeeze(sum(W(:,i,:),1)),[],1).^2> .90*sum(squeeze(sum(W(:,i,:),1)).^2,1)); % or one timepoint has >50% of the power 
end
indempty = indempty | temp;
W(:,indempty,:) = []; % Delete factors that meet the above critera

if ~isempty(H)
    H(indempty,:) = [];
end
  

end

% indi_motif = W(:,temp,:);
% for j = 1:10%size(indi_motif,2)
%    temp = squeeze(indi_motif(:,j,:));
%    temp = reshape(temp,[68 68 size(temp,2)]);
%    MotifToGif(temp,[pwd filesep sprintf('badExample%d.gif',i)],'limit',99); 
% end   