function [core_community_indices,core_community_values]= CoreCommunity(idx_louvain,Idx_knn,top_fraction,num_flag)

if nargin <3; top_fraction = 0.1; end
if nargin <4; num_flag = 0; end %use number instead of fraction

%find the motifs with the most nearest neighbors in the same graph
relationships = NaN(size(Idx_knn));
for i = 1:size(Idx_knn,1)
    for j = 1:size(Idx_knn,2)
        relationships(i,j) = idx_louvain(Idx_knn(i,j));
    end
end

%loop through each row and find those that best correspond to their
%community
number_within = NaN(size(relationships,1),1);
for i = 1:size(relationships,1)
    number_within(i) = sum(relationships(i,:)==relationships(i,1));
end

%Get only the ones that are most within community
val = unique(idx_louvain);
core_community_indices = cell(1,numel(val));
core_community_values = cell(1,numel(val));
for i = 1:numel(val)
    cur_com_id = find(idx_louvain==val(i));
    [cur_com, location] = sort(number_within(cur_com_id),'descend');  %Camden you changed this from ascend to descend on 8/10/20
    if num_flag ==1
        x = min(top_fraction,length(cur_com));
    else
        x = ceil(top_fraction*length(cur_com));
    end
    location = location(1:x);  %Camden you changed this from location(end-x+1:end) to (1:x) on 8/10/20
    core_community_indices{i} = cur_com_id(location);
    core_community_values{i} = cur_com(end-x+1:end);
end
end %function
  
 
    
    
    
    
    
    
    