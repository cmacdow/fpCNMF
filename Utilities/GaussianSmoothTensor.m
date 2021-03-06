function X = GaussianSmoothTensor(X,kernel,xydimensions,nanpxs,nobleed, verbose)
%Camden MacDowell - timeless. X is a N x K x L tensor.
%disclaimer... this is an ugly function because I'm trying to get it to work with both
%legecy and current data pipelines

if nargin <3; xydimensions = 1; end
if nargin <5; nobleed = 1; end
if nargin <6; verbose = 0; end
[N, K, L] = size(X); 


if numel(kernel) == 1 %Smooth each N across L (per K). 
    if verbose; fprintf('\nSmoothing tensor in time....');  end
    for k = 1:K
        X(:,k,:)= smoothdata(squeeze(X(:,k,:)),2,'gaussian',kernel);
    end
elseif numel(kernel) ==2 && xydimensions(1)>1 && xydimensions(2)>1
    %for imaging data. For each K converts N x L into P x P x L and spatially smooths
    if verbose; fprintf('\nSmoothing tensor in space....'); end
    for k = 1:K        
        temp = conditionDffMat(squeeze(X(:,k,:))',nanpxs,[],[xydimensions,size(X,3)]);
        temp = SpatialGaussian(temp,kernel);
        temp = conditionDffMat(temp)';
        X(:,k,:) = reshape(temp,N,1,L);        
    end
elseif numel(kernel) == 3 && xydimensions(1)>1 && xydimensions(2)>1 %for imaging data. For each K converts N x L into P x P x L and smoothes across space and time. 
    if verbose; fprintf('\nSmoothing tensor in space and time....'); end 
    X_temp = zeros(xydimensions(1)*xydimensions(2),K,L);
    for k = 1:K        
        temp = conditionDffMat(squeeze(X(:,k,:))',nanpxs,[],[xydimensions,size(X,3)]);
        bad_pxls = isnan(temp); 
        temp(bad_pxls)=0;
        temp = imgaussfilt3(temp,kernel);
        if nobleed %prevents smoothing to bleed into any masked regions (vs. smoothing over small imperfections)
            temp(bad_pxls)=0;
        end
        temp = conditionDffMat(temp)';
        X_temp(:,k,:) = reshape(temp,size(X_temp,1),1,L);
    end    
    X_temp((nanvar(reshape(X_temp,[size(X_temp,1),size(X_temp,2)*size(X_temp,3)]),[],2)<=eps)==1,:,:) = [];
    X=X_temp;    
elseif ismember(numel(kernel),[2,3]) && isempty(xydimensions)
    error('P x P dimesions for 3D conversion not supplied');
else
    error('unrecognized smoothing kernel size');
end
    
end