function data_smooth = SpatialGaussian(X,kernel,type)
%Camden MacDowell - timeless
%X = pixel x pixel x time tensor
%Applies spatial gaussian filter to all data

if nargin <3; type  = 'sigma'; end
%get nan values; 
idx = isnan(X);

%if nan values have already been removed, get the pixels with no varainces
%(e.g. nan)
if sum(idx(:))==0
    temp = reshape(X,[size(X,1)*size(X,2),size(X,3)]);
    nan_rows = nanvar(temp,[],2)<=eps;
end

%remove nanvalues for smoothing
X(isnan(X))=0; %remove nan so don't spread while smoothing

data_smooth = NaN(size(X));
for i = 1:size(X,3)
    switch type
        case 'sigma'             
            data_smooth(:,:,i) = imgaussfilt(X(:,:,i),kernel,'filterdomain','spatial');
        case 'kernel'
            data_smooth(:,:,i) = imgaussfilt(X(:,:,i),'filterdomain','spatial','filtersize',kernel);
        otherwise
            error('unknown type');
    end
    
end

data_smooth(idx)=0; %set the NaN's back to zero since they now have non-zero values

if sum(idx(:))==0
    temp = reshape(X,[size(X,1)*size(X,2),size(X,3)]);
    temp(nan_rows,:) = 0;
    data_smooth = reshape(temp,[size(X,1),size(X,2),size(X,3)]);
end

end
    