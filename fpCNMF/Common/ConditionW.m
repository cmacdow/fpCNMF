function w_out = ConditionW(w,out_dim)
%Camden MacDowell - timeless
%If out_dim is 3D then converts x to a N x K x L tensor by out_dim dims -
%   expects x to be a N x K*L 2D matrix
%If out_dim is 3D then converts x to a N*L x K matrix by out_dim dims

if numel(out_dim)==3 %convert to tensor
    K = out_dim(2); L = out_dim(3);
    w_out = zeros(out_dim);
    for lag = 0:L-1
       temp = w(:,1+K*lag:K*(lag+1));
       w_out(:,:,lag+1) = temp;
    end       
else %convert to matrix
    P = out_dim(1); K = out_dim(2); %P = N*L
    w_out = reshape(permute(w,[1,3,2]),P,K);
end

end %function
    