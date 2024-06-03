%TDT based code to cross-classify faces vs name during preparation, across
%Attention and Expectation Blocks

clear
sub_list = [1:52];
sub_out = [1,2,11,14,17,34,47,51];
sub_real = setdiff(sub_list,sub_out);
block = {'Att', 'Exp'};
levels   = {'Face','Name'};
cue_n = {'1', '2'};
count = 1;
for sub = sub_real
    for event = {'cuejit'}

              home ='E:\AttExp_fMRI\Data\derivatives\';
              sub_id = num2str(sub,'%.3d');
              subdir=['sub-' sub_id];

            resultsdir = [home subdir '/MVCC/Cue_BlockCross_Faces_vs_Name/'];
%             outdir = ['\decoding\x', model_name, dec_name, 'searchlight\', (strcat(condition(1), '_', condition(2)))];
            % Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
            cfg.results.dir = resultsdir;
            
            % Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
            beta_loc =[home subdir '/GLM_models/decoding/'];;
            
            % Set the filename of your brain mask
            cfg.files.mask = {[beta_loc '/mask.nii']};
            
            % Set the label names to the regressor names which you want to use for
            % decoding, e.g. 'button left' and 'button right'
            % don't remember the names? -> run display_regressor_names(beta_loc)

            regressor_names = design_from_spm(beta_loc);
            % since TDT only takes the first row with base names, and we
            % want to indicate block number, we have to reorder
            % regressor_names
%             regressor_names = regressor_names ([3 2 1],:);
%             if strcmp (regressor_names{1,1}(1:16), 'Sn(1) cuejit_Att')
%                 runs = {{'Sn(3)', 'Sn(5)'},{'Sn(6)', 'Sn(4)'}}
%             else
%                 runs = {{'Sn(4)', 'Sn(6)'},{'Sn(5)', 'Sn(3)'}}
%             end
                
            labelname1classA = ['*' event{1} '*' block{1,1} '*' levels{1} '*'];
            labelname1classB = ['*' event{1} '*' block{1,1} '*' levels{2} '*'];
            labelname2classA = ['*' event{1} '*' block{1,2} '*' levels{1} '*'];
            labelname2classB = ['*' event{1} '*' block{1,2} '*' levels{2} '*'];
            
            cfg = decoding_describe_data(cfg,{labelname1classA labelname1classB labelname2classA labelname2classB}, [1 2 1 2],regressor_names, beta_loc, [1 1 -1 -1]); 
            cfg.files.twoway = 1; %train and test in the two directions
            
            % cfg.design
            cfg.design = make_design_xclass(cfg); % for crossvalidation add _cv
            cfg.design.unbalanced_data = 'ok';    
            %% Set additional parameters
            cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q';
            cfg.results.overwrite = 1;
            cfg.results.output{1} = 'balanced_accuracy_minus_chance';
            cfg.results.output{2} = 'accuracy_minus_chance';
            
            cfg.plot_selected_voxels = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...
            cfg.plot_design = 1;
%     
            % Run decoding
            
            [results, cfg] = decoding(cfg);
    end
end

%% Second Level
%% A. Preprocess the resulting images (smooth and mask)

%% FLAGS
%{
The script is partly modular. You can switch specific preprocessing step on (=1) or off (=0) here.
The order of the steps is fixed at the moment though. It is assumed that you run your preprocessing
in the order: deface -> fieldmap -> realign and unwarp -> slice time correction -> coregistration -> normalization -> smoothing
%}

clear all;
do_norm=1;              % normalize images to MNI space
do_smooth=1;            % smooth images
do_mask = 1;
%% PARAMETERS
%{
specify the paramters of your dataset here
%}
% smoothing
param.smooth.fwhm = [8 8 8]; % define the smoothing kernel in mm

%% DIRECTORIES
% specify the directory in which your BIDS compatible raw data are /localized
basedir = 'E:\AttExp_fMRI\Data\derivatives';
% specify the directory to which all derivatives of the raw data are to be saved
% side note: BIDS splits your data into the raw unprocessed images, and derivatives. The latter are all images derived from the raw data, e.g. preprocessed images.
derivativesdir = 'E:\AttExp_fMRI\Data\derivatives';
% specify the directory with the results images that need to be processed
resultsdir     = {'your_results_directory'};
sub_list = [1:52];
sub_out = [1,2,11,14, 17,34,47, 51];
sub_real = setdiff(sub_list,sub_out);

prefix = 'res_balanced_accuracy_minus_chance';
%prefix = 'res_balanced_accuracy_minus_chance';
for sub = sub_real
    % convert subject number into a string for eaNsier handling
    substr=num2str(sub,'%.3d');
    subdir=['/sub-' substr];
    %% NORMALIZATION
    %normalize the images to MNI space
    if do_norm
        clear matlabbatch;
        % new normaization, using deformation fields from segmentation
        for r = 1:length(resultsdir)
            img_files{r,1} = spm_select('FPList', ...
                [derivativesdir subdir resultsdir{r}], ['^' prefix '.*\.nii$']);
        end
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList', [derivativesdir subdir '/anat/'],'^y.*.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(img_files);
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        spm_jobman('run', matlabbatch);
        clear img_files;
    end
    
    %% SMOOTHING
    % smooth the images
    if do_smooth
        clear matlabbatch;
        for r = 1:length(resultsdir)
            img_files{r,1} = spm_select('FPList', ...
                [derivativesdir subdir resultsdir{r}], ['^w' prefix '.*.nii$']);
        end
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(img_files);
        matlabbatch{1}.spm.spatial.smooth.fwhm = param.smooth.fwhm;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run', matlabbatch);
        clear matlabbatch img_files;
    end
    %% Mask image
    if do_mask == 1
        maskdir = '/GLM_models/univariate/';
        mask_file = cellstr(spm_select('FPList', ...
            [derivativesdir subdir maskdir], '^mask.nii$'));
        for r = 1:length(resultsdir)
            clear img_file img_name;
            clear matlabbatch;
            img_file = spm_select('FPList', ...
                [derivativesdir subdir resultsdir{r}], ['^sw' prefix '.*.nii$']);
            [~, img_name, ~] = fileparts(img_file);
            matlabbatch{1}.spm.util.imcalc.input = [mask_file; cellstr(img_file)];                
            matlabbatch{1}.spm.util.imcalc.output = ['masked_' img_name  '.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir = {[derivativesdir subdir resultsdir{r}]};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run', matlabbatch);
        end
    end
end
%% B. Apply Second Level using SPM12
clear matlabbatch;
cell_scans = cell(length(sub_real),1);
for s = sub_real
    cell_scans{s,1} = ['E:\AttExp_fMRI\Data\derivatives\sub-' num2str(sub_list(s),'%.3d') '\MVCC\Cue_BlockCross_Faces_vs_Name\masked_swres_balanced_accuracy_minus_chance.nii,1'];
end
cell_scans = cell_scans(~cellfun('isempty',cell_scans));
matlabbatch{1}.spm.stats.factorial_design.dir = {'E:\AttExp_fMRI\Data\derivatives\second_level\multivariate\MVCC\Cue_vs_Target/Att\Train_Cue'};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cell_scans;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'AttvExp';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
spm_jobman('run', matlabbatch);

% To plot the figure follow the following steps:
% 1. Go to SPM 
% 2. Open the results of this analysis
% 3. Select a threshold of p<0.001
% 4. Take the cluster size of the smallest significant cluster (p<0.05, FWE
% corrected)
% 5. Repeat step 3, but in cluster size use the number from step 4.
% 6. Save the result as a thresholded image.
% 7. Apply thresholded image as an overlay on MRIcroGL or similar