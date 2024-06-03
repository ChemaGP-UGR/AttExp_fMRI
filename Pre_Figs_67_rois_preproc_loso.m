%Here we extract rois based on previous analyses. To avoid double dipping,
%we apply Leave one subject out (LOSO). 
%To apply LOSO , we have to calculate the second level
% analyses including all subjects but the current one. 
%José M.G. Peñalver (cgpenalver@ugr.es)
%adapted from Ana P. Palenciano (palencianoap@ugr.es)
clear
%% Step 1: Second level analyses
% if univariate:
%spmdir = '\GLM_models\univariate\';
%contrast = con_001.nii
%if multivariate;
spmdir = '\MVPA\Cue_Att\';
contrast = 'masked_swres_balanced_accuracy_minus_chance.nii'
sub_list = [1:52];
sub_out = [1,2,11,14,17,34,47,51];
sub_real_base = setdiff(sub_list,sub_out);
do_second_level = 1; %basic step, perfom second level
do_conjunction = 1; %if the ROI of interest is the conjunction between two others
do_substraction = 0; % if the ROI of interest is the difference between two others
for sub = sub_real_base
    rois_dir = ['E:\AttExp_fMRI\Data\rois\Empirical'];
    %rois_dir2 = ['E:\AttExp_fMRI\Data\rois\Empirical']; %in case you are
    %going to need a second one
    
    rois_list = 'CUE_ATT_mask';
    %rois_list2 = 'CUE_ATT_vs_EXP_OCC_mask'; %if necessary, for conjunction calculations
    if do_second_level
        sub_real = setdiff(sub_real_base, sub); %%leave out the current subject
        %sub_real = sub_real_base;
        clear matlabbatch;
        cell_scans = cell(length(sub_real),1);
        for s = 1:length(sub_real)
            cell_scans{s,1} = ['E:\AttExp_fMRI\Data\derivatives\sub-' num2str(sub_real(s),'%.3d') spmdir contrast];
        end
        matlabbatch{1}.spm.stats.factorial_design.dir = {rois_dir};
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
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Category';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0; %here we're done with second level analyses
        matlabbatch{4}.spm.stats.results.spmmat =cellstr([rois_dir '\SPM.mat']); %show results and save them, so that we can later extract cluster size. 
        matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{4}.spm.stats.results.conspec.thresh = 0.001;
        matlabbatch{4}.spm.stats.results.conspec.extent = 0;
        matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{4}.spm.stats.results.units = 1;
        matlabbatch{4}.spm.stats.results.export{1}.ps = 1; %this will be read later to get some data
        spm_jobman('run', matlabbatch); 
        report = dir('*.ps');
        fin = fopen(report.name,'r');
        while ~feof(fin)
            s = fgetl(fin);
            if ~isempty(strfind(s,'FWEc:')) %we look for cluster size info
                index1 = strfind(s,'FWEc:') + 6;
                comas = strfind(s,','); index2 = comas(comas>index1)-1;
                FWEcluster_Threshold = s(index1:index2); %and store it here
            end
        end
        clear matlabbatch
        matlabbatch{1}.spm.stats.results.spmmat =cellstr([rois_dir '\SPM.mat']); %now we plot the results again, but applying the corresponding threshold
        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
        matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
        matlabbatch{1}.spm.stats.results.conspec.extent = str2double(FWEcluster_Threshold); %if there are several clusters that we want to split, specify here by adding a number
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.export{1}.binary.basename = rois_list; 
        spm_jobman('run', matlabbatch); 
    end

    if do_conjunction
        if exist ('matlabbatch', 'var')
            clear matlabbatch
        end
        matlabbatch{1}.spm.util.imcalc.input = {[rois_dir '\spmT_0001_' rois_list '.nii'];
                                                [rois_dir2 '\' rois_list2 '.nii']
                                                };
        matlabbatch{1}.spm.util.imcalc.output = ''; %name of the mask
        matlabbatch{1}.spm.util.imcalc.outdir = {['']}; % directory to save it
        if ~isdir(matlabbatch{1}.spm.util.imcalc.outdir{:})
            mkdir(matlabbatch{1}.spm.util.imcalc.outdir{:})
        end
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2'; %we calculate the conjunction by multiplying the two masks
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch); 
        rois_dir =  matlabbatch{1}.spm.util.imcalc.outdir{:};
        rois_list = matlabbatch{1}.spm.util.imcalc.output;
    end
    
    if do_substraction
        if exist ('matlabbatch', 'var')
            clear matlabbatch
        end
        matlabbatch{1}.spm.util.imcalc.input = {[rois_dir '\spmT_0001_' rois_list '.nii'];
                                                [rois_dir2 '\spmT_0001_' rois_list2 '.nii']
                                                };
        matlabbatch{1}.spm.util.imcalc.output = '';
        matlabbatch{1}.spm.util.imcalc.outdir = {rois_dir};
        if ~isdir(matlabbatch{1}.spm.util.imcalc.outdir{:})
            mkdir(matlabbatch{1}.spm.util.imcalc.outdir{:})
        end
        matlabbatch{1}.spm.util.imcalc.expression = 'i1-i2';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch); 
        rois_dir =  matlabbatch{1}.spm.util.imcalc.outdir{:};
        rois_list = matlabbatch{1}.spm.util.imcalc.output;
    end


home = 'E:\AttExp_fMRI\Data\';

%% Step 2: ROi preprocessing
%We select the roi for this subject and adapt it to its native space
do_correg     = 1;
do_inversNorm = 1;
do_binary     = 1;
    %% Create sbj directory
    % convert subject number into a string for easier handling
    substr = num2str(sub,'%.3d');
    substr = ['sub-' substr];
    subdir = fullfile(home,'derivatives', substr, 'anat');
    subdir2 = fullfile(home,'derivatives',substr,'func');
    roidir = fullfile(subdir,'rois');
    if ~isdir(roidir)
        mkdir(roidir);
    end
    %% Create a copy of the ROI file in the sbj folder:
    if do_conjunction || do_substraction
        roi_im = [rois_list '.nii'];
    else
        roi_im = ['spmT_0001_' rois_list '.nii'];
    end
    clear matlabbatch
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.files = ...
        {[rois_dir '\' roi_im ]};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = {roidir};
    spm_jobman('run', matlabbatch);       
    %% Inverse normalization
    if do_inversNorm
        clear matlabbatch
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = ...
            cellstr(spm_select('FPList', subdir, '^iy_rsub.*.nii'));
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample(1) = ...
            {[roidir filesep roi_im]};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
         spm_jobman('run', matlabbatch);

        %%
        % Delete the unnormlaized roi
        clear matlabbatch
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.files(1) = ...
            {[roidir filesep roi_im]};
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
        spm_jobman('run', matlabbatch);
    end

    %% Corregister (RESLICE) the ROI to the mean image (voxel size: 2.5mm)
    if do_correg
        clear matlabbatch
        matlabbatch{1}.spm.spatial.coreg.write.ref = ...
            cellstr(spm_select('FPList', subdir2, '^mean.*.nii'));
        matlabbatch{1}.spm.spatial.coreg.write.source = ...
            {[roidir filesep 'w' roi_im]};
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        spm_jobman('run', matlabbatch);

        %
        % Delete the uncoregistered roi
        clear matlabbatch
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.files(1) = ...
            {[roidir filesep 'w' roi_im]};
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
        spm_jobman('run', matlabbatch);

    end
    %% Binarize mask
    if do_binary
        clear matlabbatch;
        matlabbatch{1}.spm.util.imcalc.input  = {[roidir filesep 'rw' roi_im]};
        matlabbatch{1}.spm.util.imcalc.output = ['bin_rw' roi_im];
        matlabbatch{1}.spm.util.imcalc.outdir = {roidir};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.5';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch);

        % Delete the unbinarized roi
        clear matlabbatch
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.files(1) = ...
            {[roidir filesep 'rw' roi_im]};
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.action.delete = true;
        spm_jobman('run', matlabbatch);
    end
end