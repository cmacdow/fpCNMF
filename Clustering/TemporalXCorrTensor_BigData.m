function [tcorr_mat, lag_mat, t] = TemporalXCorrTensor_BigData(X,max_shift,verbose)
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

% precompute the shifts
fprintf('\n Precomputing and Saving Shifts')
% create a temp directory. Add random for multiple instances
rng('shuffle');
outDir = [pwd filesep sprintf('tmp_%d',randperm(10000,1)), filesep];
if ~exist(outDir,'dir'); mkdir(outDir); end
for i = 1:N
    cur_padded = cat(2, zeros(P,max_shift+1), squeeze(X(:,i,:)), zeros(P,max_shift+1));     %Zero pad
    cur_shifted = zeros(P*size(cur_padded,2),numel(t));
    for shift = 1:numel(t)
        temp = circshift(cur_padded,t(shift),2);
        temp(temp<eps)=0; %increase sparseness
        cur_shifted(:,shift) = temp(:);
    end
    fname = fullfile(outDir, sprintf('data_%08d.mat', i));
    save(fname, 'cur_shifted');
end
fprintf('.... Done')
file_list = GrabFiles('.mat',0,{outDir});

%preallocate the motifs
Z_padded = (2*(max_shift+1))+Z;
motifs = NaN(P*Z_padded,N);
for i = 1:N
   cur_padded = cat(2, zeros(P,max_shift+1), squeeze(X(:,i,:)), zeros(P,max_shift+1));     %Zero pad    
   motifs(:,i) = cur_padded(:);
end 

%compute the cross correlation
for i = 1:N %loop through each motif
    if mod(i,floor(N*.01))==0 && verbose
        fprintf('\n %d %% Complete',round(i/N*100));
    end
    
    cur_shifted = load(file_list{i},'cur_shifted');
    cur_shifted = cur_shifted.cur_shifted;
    
    %Correlate
    rho = corr(cur_shifted,motifs); 
    
    %find the best cross correlation 
    [tcorr_mat(i,:),lag_mat(i,:)] = max(rho);    
end

%delete tmp
cellfun(@(x) delete(x),file_list)
rmdir(outDir)

if verbose; fprintf('\nCalculating cross correlations took %.2g hours\n\n\n\n',toc/(60*60)); end

end