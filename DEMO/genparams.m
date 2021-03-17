classdef genparams  
    properties
        %Path Options
        local_bucket = 'Z:\';
        spock_bucket = '\jukebox\buschman\';
        repo_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\';        
        dynamic_script_path = 'Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\Spock\DynamicScripts\';
        processing_intermediates = 'Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\Scratch_Processing_Intermediates\'; %location of the intermediate files in the processing pipeline
        figure_save_directory = 'Projects\Cortical Dynamics\Parietal Cortex and Cortex Wide Dynamics\ProcessingFigures\'; %directory to save off figures during processing
        
        %dynamic script options
        sbatch_time = 59;
        sbatch_exclude = 'redshirt-n[12-49]';
        sbatch_memory = 16;
        sbatch_matlabversion = 'R2018a';
        sbatch_path = "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/";
        sbatch_name = [];               
        stack_suffix = '.ome_stack.mat' %the suffix used to take individual dff and ID for combining e.g. '_dff.mat' for corrected recs and '_dff_uncorrected.mat' for uncorrected recs
        delete_singlefiles = 1; %flag to delete individual recording sections after done concatenating
        
        %widefield post-processing parameters
        fps = 13;
        denoise_powerfrac = 0.5;
        denoise_pxlfrac = 0.01; %fraction of pixels to use. 1% seems good. Can gut check with ImpactOfPCADenoise function in gutcheck directory
        w_deconvolution = 'lucric'; %type of deconvolution (if none, then will filter by below parameters)
        w_nstd = 1; %number of standard deviations for thresholding if using w_deconvolution = 'simple_filter';
        w_filter_freq = [0.1 4]; %frequency range to filter widefield data if using no deconvolution
        w_filter_type = 'lowpass'; 
        w_normalization_method = 'full'; %pixelwise, full, or bounded
        w_norm_val = 100; %either the precentile or the value (if bounded) to normalize to
        w_chunk_dur = 120 %duration of training/testing chunks for fitting seqNMF in seconds
        w_approx_chunk_num = round(9561/(120*13)/2); %(total duration/w_chunk_dur*fps)/2 (for test and train split) This is used in pipeline to parallelize motifs fittings spock jobs without knowing the exact chunk number used. Unused will just fail as spock jobs.         
        w_pca_denoise = 0; %boolean
        
        %CNMF Defaults        
        reverse_fit = 1; 
        K = 30;
        L = 13;
        non_penalized_iter = 0;
        max_non_penalized_iter =1; 
        w_update_iter = 1;
        speed = 'fast';
        penalized_iter = 300;
        penalized_iter_refit = 100; %100
        repeat_fits = 3;
        fit_criterion = 'AIC'; %PEV, AIC, BIC;
        
        
        %specific terms for pMU
        lambda = sort(logspace(-1,-5,12), 'ascend'); %if range, then lambda sweep performed. 0.0005
        ortho_H = 1;
        ortho_W = 0;
        sparse_H = 1;
        sparse_W = 0;
      
        %General clustering Parameters
        clust_nobleed = 1; %whether to allow smoothign to bleed into masked regions
        clust_method = 'PhenoCluster';
        clust_smooth_kernel = [1,1,0.1]; %default = [1,1,0.1]
        clust_community_fraction = [0.01:0.05:1]; %if numel()>1 then will sweep and determine the best fraction per motif using the autofit_method
        autofit_method = 'autocorrelation'; %stvar or autocorrelation
        clust_removepad = 0;

        %PhenoCluster parameters
        clust_knn = 10;
        clust_louvain_restarts = 5; 

        %DBSCAN parameters
        clust_epsilon = 0.3; 
        clust_minpts = 4; 
        
        %deconvolution defaults
        d_smooth_kernel = [];
        d_rise = 200;
        d_decay = 2000; 
        d_gamma = 0.95; 
        d_smooth = 1; %default is 7
        d_kernel = 30;
        
        %miscellaneous additions
        pixel_dim = [68,68];        
        originaldimensions = [68,68];        
        verbose = 1;       
        smt_kernel = [0.1, .25, 0.5 0.75 1, 1.5]; %Use the ConfirmSmoothingLevel to check. If numel() >2 then code will choose best value. Leave blank for nosmoothing

        
    end

    methods
 
    end

end











