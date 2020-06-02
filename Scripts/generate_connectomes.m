%% GLM Abide July 12 - 2018 Sebastian Urchs - Clara Moreau
clear;

root_path = '/project/6003287/PROJECT/';
preproc_path = ['/project/6003287/PREPROCESS/'];
path_out = [root_path 'glm/'] ;
model_path = [root_path 'pheno/phenotypic_information.csv'];

files_in.networks.cambridge64 = [root_path 'template/mist_parcellation_scale_064.nii.gz'];

opt_g.min_nb_vol = 40; % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'glm_connectome'; % Specify to the grabber to prepare the files for the glm_connectome pipeline
tmp =  niak_grab_fmri_preprocess(preproc_path, opt_g);
files_in.fmri = tmp.fmri;
files_in.model.group = model_path;

opt.fdr = 0.05;
opt.folder_out = path_out; % Where to store the results


% A case-control contrast is defined here for NIAK syntax reasons in order to generate the individual seed FC matrices
% All case-control contrasts used in the paper were computed in python.
opt.test.case_control.group.contrast.Diagnosis = 1 ;
opt.test.case_control.group.contrast.Site1 = 0;
opt.test.case_control.group.contrast.Site2 = 0;
opt.test.case_control.group.contrast.FD_scrubbed = 0;
opt.test.case_control.group.flag_intercept = true;
opt.test.case_control.group.normalize_x = false;
opt.test.case_control.group.normalize_y = false;

opt.flag_test = false; 
[pipeline,opt] = niak_pipeline_glm_connectome(files_in,opt);

