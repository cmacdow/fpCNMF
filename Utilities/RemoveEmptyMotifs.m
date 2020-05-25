function [W,H] = RemoveEmptyMotifs(W,H)

if nargin <2; H =[]; end

indempty = sum(sum(W>eps,1),3)==0; % W is literally empty
indempty = indempty | (max(sum(W,3),[],1).^2> .5*sum(sum(W,3).^2,1)); % or one pixel has >50% of the power
W(:,indempty,:) = []; % Delete factors that meet the above critera

if ~isempty(H)
    H(indempty,:) = [];
end
  

end