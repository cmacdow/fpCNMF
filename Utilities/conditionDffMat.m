function [dffmat,badcols,badrows] = conditionDffMat(dff,badcols,badrows,imsize3) 

% [dffmat,badcols,badrows] = conditionDffMat(dff,badcols,badrows,imsize3) 
% converts to 3D (image format) from time x pxls format or vice-versa
% accoridng to input dff dimensions. badcols and barows are list of indices
% of cols and rows removed from time x pxls format (e.g. nan or zero from vasculature
% and mask), in case npxl ~= nX*nY
% imsize3 is a 3-element vector specifying the dimensions of the desired
% output image in case of 2d to 3d conversion
%
% LP sep 2016 - modified by CM in 2020

if nargin < 2
    badcols = [];
end
if nargin < 3
    badrows = [];
end
if nargin < 4 && size(dff,3) <= 1
    npxl    = size(dff,2)+length(badcols);
    nZ      = size(dff,1)+length(badrows);
    imsize3 = [sqrt(npxl) sqrt(npxl) nZ];
end
    
if size(dff,3) > 1
    dffmat            = reshape(dff,[],size(dff,3))';
    badcols           = find(sum(isnan(dffmat))==size(dffmat,1));    
%     badcols           = find(sum(isnan(dffmat))==size(dffmat,1) | nanvar(dffmat,[],1)<=eps);
    dffmat(:,badcols) = [];
    badrows           = find(sum(isnan(dffmat),2)==size(dffmat,2));
%     badrows           = find(sum(isnan(dffmat),2)==size(dffmat,2) | nanvar(dffmat,[],2)<=eps);
    dffmat(badrows,:) = []; 
else
    nZ     = imsize3(3);
    npxl   = imsize3(1)*imsize3(2);
    temp   = nan(nZ,npxl);
    temp(setdiff(1:nZ,badrows),setdiff(1:npxl,badcols)) = dff;
%     dffmat = reshape(temp',[nX nX nZ]);
    dffmat = reshape(temp',[imsize3(1) imsize3(2) imsize3(3)]);
end