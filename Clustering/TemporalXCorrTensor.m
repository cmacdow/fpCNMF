function [tcorr_mat, lag_mat, t] = TemporalXCorrTensor(X,max_shift,verbose)
%Camden MacDowell - timeless. 
%Compute the maximum temporal cross correlation between patterns in a
%spatio-temporal tensor.
%X is a unit x pattern x time (e.g. unit can be pixels or neurons)
%compare each pattern in X with all patterns in X and determine the maximum
%temporal x corr and the best lag.  
if nargin <3; verbose =1; end
tic 
[P, N, Z] = size(X);
assert(max_shift<=Z,'Error: max cross correlation shift is longer than the tensor pattern'); %quick catch to make sure the max shift isn't too big)
tcorr_mat = nan(N,N); %matrix of all temporal xcorrs
lag_mat = nan(N,N); %matrix of all lags
t = (-max_shift:1:max_shift); %Preallocate temporal shifts; g
if verbose; fprintf('\n\n\n\n\nCalculating Cross Correlations (now would be a good time to got get some coffee/snacks)'); end
for i = 1:N %loop through each motif
    if mod(i,floor(N*.05))==0 && verbose
        fprintf('\n %d %% Complete',round(i/N*100));
    end
    
    %Preallocate shifts and linearized motif
    cur_padded = cat(2, zeros(P,Z), squeeze(X(:,i,:)), zeros(P,Z));     %Zero pad
    cur_shifted = NaN(P*size(cur_padded,2),numel(t));
    for shift = 1:numel(t)
        temp = circshift(cur_padded,t(shift),2);
        cur_shifted(:,shift) = temp(:);
    end
    
    %best temporal cross correlation with each motifs
    for j = 1:N
       temp = cat(2, zeros(P,Z), squeeze(X(:,j,:)), zeros(P,Z));     %Zero pad       
       [tcorr_mat(i,j),lag_mat(i,j)] = max(corr(cur_shifted,temp(:)));      
    end 
end

if verbose; fprintf('\nCalculating cross correlations took %.2g hours\n\n\n\n',toc/(60*60)); end

end