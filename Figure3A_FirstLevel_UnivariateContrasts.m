% SPM12 code to estimate first level beta values on all subjects and
% perfom univariate contrasts

clear all
home = 'E:\AttExp_fMRI\Data\derivatives\'; %your derivatives folder
nRuns = [1:8]; 
sub_list = [1:52];
sub_out = [1,2,11,14,17,34, 47,51]; % subjects not included in the analyses for excesive noise, or blocks missing that prevent orthogonalization of cue pairs
sub_real = setdiff(sub_list,sub_out)

% FIRST LEVEL
for sub = sub_real
    %% Create sbj directory
    % convert subject number into a string for easier handling
    sub_id = num2str(sub,'%.3d');
    subdir=['sub-' sub_id];
    
    for analysis = {'univariate'}
        
        destination = [home subdir '\GLM_models\' analysis{1} '\subsampled'];
        
        if ~isdir(destination)
            mkdir(destination);
        end
        
        if strcmp(analysis{1},'univariate') %data undergoes different preprocessing steps depending if betas are estimated for decoding or univariate analyses
            prefix_im   = 'swau';
        elseif strcmp(analysis{1},'decoding')
            prefix_im   = 'au';
        end
        
        %-----------------------------------------------------------------------
        % Job saved on 03-Dec-2020 15:46:21 by cfg_util (rev $Rev: 7345 $)
        % spm SPM - SPM12 (7487)
        % cfg_basicio BasicIO - Unknown
        %-----------------------------------------------------------------------
        clear matlabbatch;
        matlabbatch{1}.spm.stats.fmri_spec.dir = {destination};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.73;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25;        
        %%
        nRuns = [1:8]; %number of runs

        for r = nRuns
          matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans =  ...
                cellstr(spm_select('FPList',[home '\sub-' sub_id '\func\'],...
                ['^' prefix_im 'sub.*run-' num2str(runsIn(r)) '_bold.nii$']));
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}); % you have to fill these fields with your own values
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = ...
                {[home 'sub-' sub_id '\func\sub-' sub_id '_task-attexp_run-' num2str(runsIn(r)) '_events_subsampled.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = ...
                {[home 'sub-' sub_id '\func\rp_sub-' sub_id '_task-attexp_run-' num2str(runsIn(r)) '_bold.txt']}; %regresores de movimiento
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        end
        
        %%
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        spm_jobman('run', matlabbatch);
        
        %%
        clear matlabbatch;
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[ destination '\SPM.mat']};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', matlabbatch);
        
        %%
        disp([sub_id ': Response-related ' analysis{1} ' model estimated!']);
    end 
end

%CONTRASTS
for sub = sub_real
        % convert subject number into a string for easier handling
        sub_id = num2str(sub,'%.3d');
        subdir=['sub-' sub_id];

        glm_dir = [home subdir '/GLM_models/univariate/'];
        load([glm_dir filesep 'SPM.mat']);

        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = {[glm_dir '/SPM.mat']};

        %% Figure 3A
        % Cue Att - Cue Exp 
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Cue Att - Cue Exp';
        wv = zeros(size(SPM.xX.name));
        wv(contains(SPM.xX.name,'Att') & contains(SPM.xX.name,'cuejit')) = 1;
        wv(contains(SPM.xX.name,'Exp') & contains(SPM.xX.name,'cuejit')) = -1;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = wv;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        
        % Cue Exp - Cue Att 
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Cue Exp - Cue Att';
        wv = zeros(size(SPM.xX.name));
        wv(contains(SPM.xX.name,'Exp') & contains(SPM.xX.name,'cuejit')) = 1;
        wv(contains(SPM.xX.name,'Att') & contains(SPM.xX.name,'cuejit')) = -1;
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = wv;
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    

        matlabbatch{1}.spm.stats.con.delete = 1;
        spm_jobman('run', matlabbatch); 

end

%% Second Level

%Cue Att - Cue Exp
 for s = 1:length(sub_real)
            cell_scans{s,1} = ['E:\AttExp_fMRI\Data\derivatives\sub-' num2str(sub_real(s),'%.3d') '\GLM_models\univariate\con_0001.nii,1'];
 end
    matlabbatch{1}.spm.stats.factorial_design.dir = {['E:\AttExp_fMRI\Data\derivatives\second_level\' type '\' confold{c} filesep conname{c} ]};
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
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Cue Att - Cue Exp';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;
    spm_jobman('run', matlabbatch);
    
 % Cue Exp - Cue Att
 
 for s = 1:length(sub_real)
            cell_scans{s,1} = ['E:\AttExp_fMRI\Data\derivatives\sub-' num2str(sub_real(s),'%.3d') '\GLM_models\univariate\con_0002.nii,1'];
 end
    matlabbatch{1}.spm.stats.factorial_design.dir = {['E:\AttExp_fMRI\Data\derivatives\second_level\' type '\' confold{c} filesep conname{c} ]};
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
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Cue Exp - Cue Att';
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