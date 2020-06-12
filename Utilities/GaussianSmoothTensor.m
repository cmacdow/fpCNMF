function X = GaussianSmoothTensor(X,kernel,xydimensions,nanpxs)
%Camden MacDowell - timeless. X is a N x K x L tensor. 

if nargin <3; xydimensions = 1; end

[N, K, L] = size(X); 


if numel(kernel) == 1 %Smooth each N across L (per K). 
    fprintf('\nSmoothing tensor in time....');    
    for k = 1:K
        X(:,k,:)= smoothdata(squeeze(X(:,k,:)),2,'gaussian',kernel);
    end
elseif numel(kernel) ==2 && xydimensions(1)>1 && xydimensions(2)>1
    %for imaging data. For each K converts N x L into P x P x L and spatially smooths
    fprintf('\nSmoothing tensor in space....');
    for k = 1:K        
        temp = conditionDffMat(squeeze(X(:,k,:))',nanpxs,[],[xydimensions,size(X,3)]);
        temp = SpatialGaussian(temp,kernel);
        temp = conditionDffMat(temp)';
        X(:,k,:) = reshape(temp,N,1,L);        
    end
elseif numel(kernel) == 3 && xydimensions(1)>1 && xydimensions(2)>1
    %for imaging data. For each K converts N x L into P x P x L and smoothes across space and time. 
    fprintf('\nSmoothing tensor in space and time....');
    for k = 1:K        
        temp = conditionDffMat(squeeze(X(:,k,:))',nanpxs,[],[xydimensions,size(X,3)]);
        bad_pxls = isnan(temp); %assumes images are masked. Saves this infromation so that smoothing doesn't bleed into masked regions
        temp(bad_pxls)=0;
        temp = imgaussfilt3(temp,'FilterSize',kernel);
        temp(bad_pxls)=NaN;
        temp = conditionDffMat(temp)';
        X(:,k,:) = reshape(temp,N,1,L);
    end 
elseif ismember(numel(kernel),[2,3]) && isempty(xydimensions)
    error('P x P dimesions for 3D conversion not supplied');
else
    error('unrecognized smoothing kernel size');
end
    
end