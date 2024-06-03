%Code to apply voxel selectivity ranking, plotted in Figure 6. 
%José M.G. Peñalver (cgpenalver@ugr.es)
%adapted from Carlos González García (cgonzalez@ugr.es)

clear

home = 'E:\AttExp_fMRI\Data\derivatives\';
cd(home);
nRun = 8;
sub_list = [1:52];
sub_out = [1,2,11,14,17,34,47,51];
sub_real = setdiff(sub_list,sub_out);
block = {'Att', 'Exp'};
validity = {'Val','Inv'};
roi = {'spmT_0001_CAtt_FN_mask'; 'spmT_0001_CExp_FN_mask';'CTmvcc_Att_FN_mask_VVC'; 'CTmvcc_Att_FN_mask_FO'; 'CTmvcc_Exp_FN_mask' }; %these are the ROIs obtained from the previous LOSO (pre_Figs_67_rois_preproc_loso.m)

for r = 1:length(roi) 
        sub_count= 1;
    for sub = sub_real
    for bl = 1:2
    for event = {'cue'}
        levels   = {'Face','Name'};
            %% Create sbj directory
            
            sub_id = num2str(sub,'%.3d');
            subdir=['sub-' sub_id];
            
            beta_loc = [home subdir '/GLM_models/decoding/subsampled'];
            %% Initialize cfg variable:
            cfg = decoding_defaults;
           
            % Mask to restrict the analysis to the participant's brain
            % (native space, unsmoothed)
            cfg.files.mask = [home subdir '\anat\rois\' roi{r} '.nii'];
            
            % The following function extracts all beta names and corresponding run
            % numbers from the SPM.mat
            regressor_names = design_from_spm(beta_loc);
            labelname1 = [event{1} '*' block{1,bl} '*' levels{1} '*'];
            labelname2 = [event{1} '*' block{1,bl} '*' levels{2} '*'];
            % Extract all information for the cfg.files structure (labels will be [1 -1] )
            cfg = decoding_describe_data(cfg,{labelname1 labelname2},[1 -1],regressor_names,beta_loc);
            % load voxel info per condition
            data = decoding_load_data(cfg);
            res (1,:) = data.data(1,:);
            res (2,:) = data.data(2,:);
            res (3,:) = data.data(3,:);
            res (4,:) = data.data(4,:);
            res (5,:) = data.data(5,:);
            res (6,:) = data.data(6,:);
            res (7,:) = data.data(7,:);
            res (8,:) = data.data(8,:);          
  
        res = res';
        for cond = 1:length (res)
            res_sorted_target(cond,:) = sort (res(cond,:));
            for c = 1:8
                res_ind(cond,c) = find(res_sorted_target(cond,c) == res(cond,:));
            end
        end
        clear res
         
    end

    for event = {'target'}
        levels   = {'Face','Word'};
            %% Create sbj directory
            
            sub_id = num2str(sub,'%.3d');
            subdir=['sub-' sub_id];
            %CHANGE before running analyses
            resultsdir = 'E:\AttExp_fMRI\Results\vox_order'
            if ~isfolder(resultsdir)
                mkdir(resultsdir);
            end
            cd(resultsdir);
            
            %% Initialize cfg variable:
            cfg = decoding_defaults;
            
            %% cfg.files: .mask; .name; .descr; .chunk; .label
            
            beta_loc = [home subdir '/GLM_models/decoding/subsampled'];
            
            % Mask to restrict the analysis to the participant's brain
            % (native space, unsmoothed)
            cfg.files.mask = [home subdir '\anat\rois\' roi{r} '.nii'];
            
            % The following function extracts all beta names and corresponding run
            % numbers from the SPM.mat
            regressor_names = design_from_spm(beta_loc);
            labelname1 = [event{1} '*' block{1,bl} '*' levels{1} '*'];
            labelname2 = [event{1} '*' block{1,bl} '*' levels{2} '*'];
            % Extract all information for the cfg.files structure (labels will be [1 -1] )
            cfg = decoding_describe_data(cfg,{labelname1 labelname2},[1 -1],regressor_names,beta_loc);
            % load voxel info per condition
            data = decoding_load_data(cfg);

            res (1,:) = mean(data.data(1:2,:));
            res (2,:) = mean(data.data(3:4,:));
            res (3,:) = mean(data.data(5:6,:));
            res (4,:) = mean(data.data(7:8,:));
            res (5,:) = mean(data.data(9:10,:));
            res (6,:) = mean(data.data(11:12,:));
            res (7,:) = mean(data.data(13:14,:));
            res (8,:) = mean(data.data(15:16,:)); 
        res = res';
        for cond = 1:length (res)
            res_sorted_target(cond,:) = res(cond,(res_ind(cond,:)));
        end
        clear res
        m_res_sorted_target = mean(res_sorted_target);
        if bl==1
        One_att (sub_count,:) = m_res_sorted_target(1);
        Two_att (sub_count,:) = m_res_sorted_target(2);
        Three_att (sub_count,:) = m_res_sorted_target(3);
        Four_att (sub_count,:) = m_res_sorted_target(4);
        Five_att (sub_count,:) = m_res_sorted_target(5);
        Six_att (sub_count,:) = m_res_sorted_target(6);
        Seven_att (sub_count,:) = m_res_sorted_target(7);
        Eight_att (sub_count,:) = m_res_sorted_target(8);
        
        else
        One_exp (sub_count,:) = m_res_sorted_target(1);
        Two_exp (sub_count,:) = m_res_sorted_target(2);
        Three_exp (sub_count,:) = m_res_sorted_target(3);
        Four_exp (sub_count,:) = m_res_sorted_target(4);
        Five_exp (sub_count,:) = m_res_sorted_target(5);
        Six_exp (sub_count,:) = m_res_sorted_target(6);
        Seven_exp (sub_count,:) = m_res_sorted_target(7);
        Eight_exp (sub_count,:) = m_res_sorted_target(8);
        end
        
        
    end
        if bl == 1
            sorted_at(:,sub_count,r) = m_res_sorted_target;
        else
            sorted_ex(:,sub_count,r) = m_res_sorted_target;
        end
    end
            % fit linear regression to the ranked parameter estimates
            % (polyfit computes a least squares polynomial for the data,
            % generates the coefficients of the polynomial, which will be then
            % used to model a straight line (polynomial degree = 1) to fit the data
            x_images = (1:size(sorted_at,1))';
            fit_at = polyfit(x_images,sorted_at(:,sub_count,r),1); % first column is the slope coefficient
            slope_at(sub_count,r) = fit_at(1); %valor del coeficiente de la pendiente
            fit_ex = polyfit(x_images,sorted_ex(:,sub_count,r),1);
            slope_ex(sub_count,r) = fit_ex(1);
            
            % evaluate the polynomial
            % (polyval evaluates the polynomial in the data. It generates a
            % straight line (or curve if polynomial degree > 1) to fit the data
            % based on the coefficients found using polyfit
            at_eval(:,sub_count,r) = polyval(fit_at,x_images);
            ex_eval(:,sub_count,r) = polyval(fit_ex,x_images);
    sub_count= sub_count+1;
    clear m_res_sorted_target
    end
   
    res_table{r} = table (One_att, Two_att, Three_att, Four_att, Five_att, Six_att, Seven_att, Eight_att,...
        One_exp, Two_exp, Three_exp, Four_exp, Five_exp, Six_exp, Seven_exp, Eight_exp);
    %res_m = [res.One, res.Two, res.Three, res.Four, res.Five, res.Six, res.Seven, res.Eight];
%     writetable(res_table{4}, [savedir '\mvcc_Exp_FN_results.csv']);
%     writetable(res_table{4}, [savedir '\mvcc_Exp_FN_results.xlsx']);

    
end
 %% statistical analysis
    %% check if assumptions are met
    % before comparing the slopes of pre and post, two checks are made in
    % Ritcher deLange JNS (2018). The first one is comparing each slope to 0
    % (one-tail one-sample ttest; same results with two-tails) to make sure that there is a positive
    % regression slope, which would mean that the ranking of the gray scale
    % images generelizes to pre/post images. I perform this here:
    slope_at(slope_at==0) = NaN; % get rid of nans
    [h_at,p_at,ci,stats_at] = ttest(slope_at,0,'Tail','right');
    slope_ex(slope_ex==0) = NaN; % get rid of nans
    [h_ex,p_ex,ci,stats_ex] = ttest(slope_ex,0,'Tail','right');


    % %version with wilkoxon
    for r = 1:length(roi)
        [p_at(1,r),h_at(1,r)] = signrank(slope_at(:,r),0,'Tail','right');
        [p_ex(1,r),h_ex(1,r)] = signrank(slope_ex(:,r),0,'Tail','right');
    end
    [at_signif,at_index_signif] = fdr(p_at,0.05); % FDR-correct for multiple comparisons
    [ex_signif,ex_index_signif] = fdr(p_ex,0.05); % FDR-correct for multiple comparisons
    % from this analysis, we can conclude that in X rois, slope coefficients
    % differ from zero
    
%% compare the slopes of pre and post
%for each ROI, perform a two-tailed paired t-test. If expectation
%suppresion (i.e. pre - post) scales with voxel selectivity (i.e.
%dampening), then the slope of pre should be higher. If the opposite is
%found, this would be in line with the sharpening account.
[h_comp,p_comp,ci,t_comp] = ttest(slope_at,slope_ex);
%wilkoxon
for roi = 1:2
[p_comp(1,roi),h_comp(1,roi)] = signrank(slope_at(:,roi),slope_ex(:,roi));
end
[signif_comp,signif_idx_comp] = fdr(p_comp(1:2),0.05); % FDR-correct for multiple comparisons

%% SAVE
savedir = 'E:\AttExp_fMRI\Results\vox_order';
save([savedir '/results_att_split'], 'signif_comp', 'signif_idx_comp', 'p_comp', 'slope_at', 'slope_ex', 'p_at', 'p_ex', 'h_at', 'h_ex', 'at_eval', 'ex_eval', 'res_table')

% to plot the results as in Figure 6, plot slope_ex or slope_at for the
% slope, and res_table for the individual condition data points