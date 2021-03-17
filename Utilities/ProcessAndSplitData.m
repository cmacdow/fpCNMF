function [data_norm, nanpxs, data_train, data_test] = ProcessAndSplitData(fn,save_fn,parameter_class)
%Camden MacDowell - timeless
%filters, normalizes, and splits data in fn into training and test test
%fn can be the full file path or a stack and opt structure. 
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

gp = loadobj(feval(parameter_class)); 

%% SUPER UGLY CONTIGENCIES DUE TO LEGACY DATA
if ischar(fn) %load data
    warning('Camden you may want to check how things are transposed')
    fprintf('\n\tLoading data')
    %load data
    temp = load(fn);
    data = temp.dff;
    opts = temp.opts;
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
        
    otherwise
        error('unknown w_deconvolution option. Check general parameters'); 
end

%Denoise with PCA (removed banded pixels)
if gp.w_pca_denoise
    data = conditionDffMat(data,nanpxs,[], [x,y,z]);
    data = DenoisePCA(data);
    [data,~] = conditionDffMat(data);
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
    save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','-v7.3')
    fprintf('\n\tDONE')
end


%% for saveing off figure use; 
% [path, name] = fileparts(save_fn);
% saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),[path filesep name, '_ProcessAndSplitFigures'],0); %close all;


end















