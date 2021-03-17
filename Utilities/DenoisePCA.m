function [stack_denoise,handles, gp] = DenoisePCA(stack,varargin)
%Camden MacDowell. Initially denoise imaging data with PCA. This is best
%used to remove single spot noise pixels but can also be used to manually
%remove artifacts. 

gp = general_params;
gp = ParseOptionalInputs(gp,varargin);

[x,y,z] = size(stack);

[stack_flat,nanpxs] = conditionDffMat(stack); 
[coef, score, ~, ~, ~, mu] = pca(stack_flat');

%find components in which x pixels has over x percent of the power
num_pxls = round(gp.denoise_pxlfrac*size(score,1));
bad_comp = sum(maxk(score,num_pxls,1).^2) > gp.denoise_powerfrac*sum(score.^2,1);

%Create figure of the pixels corresponding to those components
[~,bad_idx] = maxk(score(:,bad_comp),num_pxls,1);
mask = ones(size(score,1),1);
mask(bad_idx)=2; 
mask = conditionDffMat(mask, nanpxs,[],[x,y,1]);
figure; imagesc(mask); axis equal; title('PCA Denoised Pixels in Components that Are Removed');

%Remove those components
stack_denoise = score(:,bad_comp==0)*coef(:,bad_comp==0)'+repmat(mu,size(score,1),1);

% %make figure showing the comparison; 
% figure('position',[75,558,1787,420]); hold on; 
% subplot(1,3,1); imagesc(stack_flat',[0 0.5]); title('orig');
% subplot(1,3,2); imagesc(stack_denoise,[0 0.5]); title('denoise');
% subplot(1,3,3); imagesc(stack_flat'-stack_denoise,[-0.5 0.5]); title('orig-denoise');
% handles = get(groot, 'Children'); 

fprintf('\n\t PCA Denoising Removed %d Components',sum(bad_comp));

%This can lead to a small number of frames (<1%) to contain negative values. Set those to zero 
stack_denoise(stack_denoise<0)=0;

%recondition
stack_denoise = conditionDffMat(stack_denoise',nanpxs,[],[x,y,z]);

end
