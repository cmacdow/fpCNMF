function Y = MaskTensor(X,nanpxs,dims)
%camden macdowell - timeless

if size(X,1) < dims(1)
    Y = zeros(dims); 
    temp = ones(dims(1),1);
    temp(nanpxs)=0;
    Y(temp==1,:,:) = X;
else
    Y = X;
    Y(nanpxs,:,:) = [];
end

end
