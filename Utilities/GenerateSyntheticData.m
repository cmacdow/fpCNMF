function [data, W, fh] = GenerateSyntheticData(varargin)
%% parse optional inputs
opts.blocks = 1;
opts.N = 20;
opts.T = 300;
opts.L = 15;
opts.K = 4;
opts.frac_active_time = 0.01;
opts.no_overlap_W = 1; 
opts.no_overlap_H = 1; 
opts.sigma_W = 0;
opts.sigma_data = 0.05;
opts.verbose = 1; 

opts = ParseOptionalInputs(opts,varargin);

if opts.T<20*opts.L
    warning(' T < 20L may lead to some wacky synthetic data')
end

%% Generate synthetic data
%Generate Motifs
W = SyntheticW(opts.K,opts.N,opts.L,opts.no_overlap_W);

data = cell(1,opts.blocks);
for i = 1:opts.blocks %option to build multiple blocks of data with slightly different motifs between blocks
   %Optionally add Noise to Motifs
   noise = normrnd(0,opts.sigma_W,opts.N,opts.K,opts.L);
   noise(noise<0)=0; %prevent noise from adding negative components; 
   W_noisy = W + noise;

   %Generate Temporal Weightings
   H = SyntheticH(opts.K,opts.T,opts.L,opts.frac_active_time,opts.no_overlap_H);

   %Build full data and optionally add noise
   noise = normrnd(0,opts.sigma_data,opts.N,opts.T);
   noise(noise<0)=0; %prevent noise from adding negative components;
   data{i} = tensor_convolve(W_noisy,H)+noise;  
   
   %optional visualization of the example block
   if opts.verbose && i==1; fh=VisualizeData(data{i},W_noisy,H); end 
   
end

if opts.verbose ==0; fh = []; end

end

