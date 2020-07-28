function [idx, x_sort] = SortMotifs(x,method)
%A quick shell for the sort function to keep things consistent across data
%analysis code. Also is easy to add new sorting methods 
%X is (n x motif) list of loadings where n is trials/mice/days/etc.
if nargin <2; method = 'mean'; end

switch method
    case 'mean'
        [~,idx] = sort(nanmean(x,1),'descend');
    case 'median'
        [~,idx] = sort(nanmedian(x,1),'descend');
    otherwise 
        error('Unknown averaging method');
end

x_sort = x(:,idx);

end %function