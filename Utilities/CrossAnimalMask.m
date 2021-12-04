function [w_mask, nanpxs_all] = CrossAnimalMask(W,nanpxs)
%Camden MacDowell - timeless
%Mask all data in W using only pixels shared across all animals
%W and nanpxs are cell arrays

%universal mask
nanpxs_all = unique([nanpxs{:}]);

%reconditional all 
w_mask = cell(1,numel(W));
for i = 1:numel(W)
    for j = 1:size(W{i},2)
        temp = conditionDffMat(squeeze(W{i}(:,j,:))',nanpxs{i});
        [x,y,z] = size(temp);
        temp = reshape(temp,[x*y,z]);
        temp(nanpxs_all,:) = NaN;
        temp = reshape(temp,x,y,z);
        w_mask{i}(:,j,:) = conditionDffMat(temp)';
    end
end

end %function end