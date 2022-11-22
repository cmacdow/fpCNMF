function [data_norm, nanpxs, data_train, data_test] = ProcessAndSplitData(fn,save_fn,parameter_class)
%Camden MacDowell - timeless
%filters, normalizes, and splits data in fn into training and test test
%fn can be the full file path or a stack and opt structure. 
if ~ispc
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'))
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

fprintf('splitting data');

gp = loadobj(feval(parameter_class)); 

%% Ugly contingencies for legacy data
if ischar(fn) %load data
    warning('Camden you may want to check how things are transposed')
    fprintf('\n\tLoading data')
    %load data
    temp = load(fn);
    data = temp.dff;
    opts = temp.opts;
    %reshape if conditioned
    if size(data,3)==1
        data = conditionDffMat(data,temp.nanpxs);
    end
    clear temp;
elseif isstruct(fn)
    data = fn.dff;
    opts = fn.opts;
else
    data = fn;
    opts = gp;
end
%%

%condition data and remove nan pxls
[x,y,z] = size(data);
[data,nanpxs] = conditionDffMat(data); %nanpxs are the same on each iteration so fine to overwrite

%filter data
switch gp.w_deconvolution
    case 'filter_thresh'
        fprintf('\n\tFiltering and thresholding data')
        data = filterstack(data, opts.fps, gp.w_filter_freq, gp.w_filter_type, 1, 0);
        %Remove negative component of signal and find bursts as in MacDowell 2020. Legacy. %Deconvolution with non-negative output is preferred (maintains more data, comes with own assumptions). 
        for px = 1:size(data,2)
           temp = data(:,px);
           temp(temp<=nanmean(temp(:))+gp.w_nstd*std(temp(:))) = 0;
           data(:,px) = temp;
        end
        
        data(data<(nanmean(data(:))+gp.w_nstd*nanstd(data(:))))=0;
    case 'lucric'
        fprintf('\n\tPerforming a Lucy-Goosey Deconvolution (Lucy-Richardson)\n')
        for px = 1:size(data,2)
           data(:,px) = lucric(data(:,px),gp.d_gamma,gp.d_smooth,gp.d_kernel);
        end
    case 'only_filter'
        fprintf('\n\tFiltering data')
        data = filterstack(data, opts.fps, gp.w_filter_freq, gp.w_filter_type, 1, 0);
    case 'fNN' %feedforward neural network trained per animal
        fprintf('loading pretrained neural network');
        %load train nn(see GeneratefNN_spock.sh)        
        if ispc 
            load(['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
            temp = load('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\PreprocessedImaging\restingstate_processed_fn.mat','dff_list');
            net_indx = find(strcmp(temp.dff_list,fn)==1);
        else
            load(['/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging' filesep gp.fit_nn_fn],'trained_opts');
            temp = load('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/PreprocessedImaging/restingstate_processed_fn.mat','dff_list_bucket');
            net_indx = find(strcmp(temp.dff_list_bucket,fn)==1);            
        end       
        params = trained_opts{net_indx}.feedforwardparams;
        net =trained_opts{net_indx}.shallowfeedforward;
        %loop through each pixel
        for px = 1:size(data,2)            
            xtemp = createRollingWindow(data(:,px), params.win)'; %t-n:t-1        
            stPred = net(xtemp)';        
            %NaN pad to match timepoints
            data(:,px)=NaN;
            data(ceil(params.win/2):end-floor(params.win/2),px)=stPred;
        end
        %remove pad
        padidx = sum(isnan(data),2)==size(data,2); 
        data(padidx,:)=[];
        z =  size(data,1);
    case 'none'
        fprintf('doing nothing');
            
    otherwise
        error('unknown w_deconvolution option. Check general parameters'); 
end

%Denoise with PCA (removed banded pixels)
if gp.w_pca_denoise 
    fprintf('denoising with pca');
    data = conditionDffMat(data,nanpxs,[], [x,y,z]);
    data = DenoisePCA(data);
    [data,~] = conditionDffMat(data);
end

%optionally bin data
if gp.temporalbin
   fprintf('temporally binning')
   if isEven(z)
       data = data(1:2:end,:)+data(2:2:end,:);
   else
       data = data(1:2:end-1,:)+data(2:2:end,:);
   end
   opts.fps = floor(opts.fps/2);
   z=size(data,1);
end

%normalize to 0 to 1 
fprintf('\n\tPerforming %s normalization to %d value', gp.w_normalization_method, gp.w_norm_val);
switch gp.w_normalization_method
    case 'pixelwise' %each between x and xth pixel intensity
        data_norm = NaN(size(data));
        for px = 1:size(data,2)
            data_norm(:,px) = (data(:,px))/(prctile(data(:,px),gp.w_norm_val));
        end             
    case 'full' %normalize using the percentile of the maximum         
        data_norm = data/prctile(data(data>eps),gp.w_norm_val);          
    case 'bounded'
        data_norm = (data)/(gp.w_norm_val(2)); %normalize between zero and the upper bound     
    case 'none'
        data_norm = data;
    otherwise
        error('Unknown normalization method. Check general params')
end

%transpose (fpCNMF operates rowwise)
data_norm = data_norm';

%Chunk (since MU is suboptimal for >4500 timepoints. As many as possible
num_chunks = floor(z/(gp.w_chunk_dur*opts.fps));
if ~isEven(num_chunks)% need even number
    num_chunks = num_chunks-1; 
end

%remove the remainder and reshape into chunks
data_trim = data_norm(:,1:end-mod(z,num_chunks*(gp.w_chunk_dur*opts.fps)));

data_chunked = cell(1,num_chunks);
for i = 1:num_chunks 
    data_chunked{i} = data_trim(:,1+(i-1)*(gp.w_chunk_dur*opts.fps):(i*(gp.w_chunk_dur*opts.fps)));
end

%alternate testing and training
data_train = cat(3,data_chunked{1:2:num_chunks});
data_test = cat(3,data_chunked{2:2:num_chunks});

%save off the data in the scratch directory and the nanpxs
if ~isempty(save_fn)
    fprintf('\n\tSaving data')
    if strcmp(gp.w_deconvolution,'fNN')
        save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','padidx','-v7.3')    
    else %legacy
        save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','-v7.3')
    end
    fprintf('\n\tDONE')
end


%% for saveing off figure use; 
% [path, name] = fileparts(save_fn);
% saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),[path filesep name, '_ProcessAndSplitFigures'],0); %close all;


end















