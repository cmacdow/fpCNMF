function [train, test] = SplitData(data,dur)
%Camden MacDowell - timeless
%breaks data (N x T) into equal chunks of N x dur. Splits equally into test
%and train. Surplus time points or extra chunks are excluded. 

[~,T] = size(data);

num_chunks = floor(T/dur);
if mod(num_chunks,2)% need even number to equally split test/train
    num_chunks = num_chunks-1; 
end

%remove the remainder and reshape into chunks
data_trim = data(:,1:end-mod(T,num_chunks*dur));
data_trim = reshape(data_trim,[size(data_trim,1),num_chunks,dur]);

%alternate testing and training. 
% train = data_trim(:,1:2:num_chunks,:);
% test = data_trim(:,2:2:num_chunks,:);
train = arrayfun(@(x) squeeze(data_trim(:,x,:)), 1:2:num_chunks,'UniformOutput',0);
test = arrayfun(@(x) squeeze(data_trim(:,x,:)), 2:2:num_chunks,'UniformOutput',0);




end %function end