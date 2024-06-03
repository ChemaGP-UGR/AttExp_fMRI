% José M. G. Peñalver (cgpenalver@ugr.es)
%adapted from Hebart et al. (2018) 
% Commonality analysis using Pearson correlations
clear
addpath('helper_functions')

% In hebart's paper they average both modalities
mean_fmri = 1;
mean_meeg = 1;
nperms =  5000; 
save_perms = 1;
alpha = 0.05;

resultsdir = 'E:\AttExp_fMRI\Results\Fusion\Full\commonality\cue\emp\block_cat_cue';
if ~isdir (resultsdir)
    mkdir (resultsdir)
end
%% Load data first
% EEG
temp = load('E:\AttExp_fMRI\Results\Fusion\EEG\cue\results\LDC\rdms\result_ldc_empirical.mat');
for i = 1:size(temp.RDMs,4)
    RDMs{i} = temp.RDMs(:,:,:,i);
end
meeg_rdms = RDMs';
t = temp.cfg.tm.times;

% fMRI
roi = {'spmT_0001_CAtt_FN_mask'; 'spmT_0001_CExp_FN_mask';'CTmvcc_Att_FN_mask_FO'; 'CTmvcc_Att_FN_mask_VVC'; 'CTmvcc_Exp_FN_mask'};%'V1_bilateral'; 'SPC_bilateral'; 'M1_bilateral'; 'LOC_bilateral'; 'IFS_bilateral'; 'FG_bilateral'
for r = 1:length(roi)
temp = load(['E:\AttExp_fMRI\Results\Fusion\fMRI\cue\results\LDC\rdms_the\' roi{r} '\result_ldc.mat']);

RDM_fMRI(:,r) = temp.result';
end
ROIs = roi;


%% Set up model variables

%models are defined in the RSA functions for both EEG and fMRI. Here we
%just outline them. 

%if cue
cat_model = [repmat([0,0,0,0,1,1,1,1],4,1);repmat([1,1,1,1,0,0,0,0],4,1)];
cue_model = [repmat([0,1,0,1;1,0,1,0],2,1),ones(4);ones(4),repmat([0,1,0,1;1,0,1,0],2,1)];

%% Select lower triangular matrix
% convert to vectors
RDV_category = mysquareform(cat_model);
RDV_cue = mysquareform(cue_model);
%wait to do the same with EEG and fMRI matrices

%% Average fmri / eeg

if mean_fmri
    for j = 1:size(RDM_fMRI,1)
        mean_data = zeros(8);
        for i = 1:size(RDM_fMRI,2)
            mean_data = RDM_fMRI{j,i}+mean_data;
        end
        mean_data = mean_data/size(RDM_fMRI,2);
        pre_RDM_fMRI(:,:,j) = mean_data;
    end
end

if mean_meeg
    mean_data = zeros(8,8, 142);
    for i = 1:size(meeg_rdms,1)
        mean_data = meeg_rdms{i}+mean_data;
    end
    mean_data = mean_data/size(meeg_rdms,1);
    pre_RDM_EEG = mean_data;  
end

% smooth EEG data
for time = 1:length (pre_RDM_EEG)
    if time>5 && time<137
        pre_RDM_EEG_s(:,:,time) = mean(pre_RDM_EEG(:,:,time-2:time+2), 3);  %this applies a 5 second moving average window
    else
        pre_RDM_EEG_s(:,:,time) = pre_RDM_EEG(:,:,time);
    end
end

RDM_EEG = pre_RDM_EEG_s;
RDM_fMRI = pre_RDM_fMRI;
  
RDV_fMRI = mysquareform(RDM_fMRI);
RDV_EEG = mysquareform(RDM_EEG);

%% This would be simple MEG-fMRI fusion

Fusionmat(:,:) = corr(RDV_fMRI,RDV_EEG,'type','pearson');


%% We want more complicated model-based fusion (based on commonality analysis)

disp('Running commonality analysis...')
for j_roi = length(roi):-1:1
    
    xMRI = RDV_fMRI(:,j_roi);
    xcat = RDV_category;
    xcue = RDV_cue;
    
    for i_time =1:size(RDM_EEG,3)

        y = RDV_EEG(:,i_time);
        
        % and now again with another calculation

        rEEG_MRIcuecat= correlate([y xMRI xcue xcat],'type','pearson','method','semipartialcorr');
        
        rEEG_MRIcat = correlate([y xMRI xcat],'type','pearson','method','semipartialcorr');
        rEEG_MRIcue = correlate([y xMRI xcue], 'type','pearson','method','semipartialcorr');
        
        CMEGMRIcat(i_time,j_roi) = abs(rEEG_MRIcue(2,1).^2-rEEG_MRIcuecat(2,1).^2);
        CMEGMRIcue(i_time,j_roi) = abs(rEEG_MRIcat(2,1).^2-rEEG_MRIcuecat(2,1).^2); 
    end
end
disp('done.')

%% Now we run a randomization test, randomizing the rows and columns of the MEG data

% Since all values in MEG are unique, we don't need to use tiedrank
% and since all other values don't change in the 612 steps and the 1000
% permutations, we can just calculate them once, rank transform them, and
% instead of running spearman correlation run pearson correlation.

% if you want to speed this up, run variance partitioning on multiple
% permutations in parallel (treating them as separate), on the ranks of the
% permutations

if exist([resultsdir '\permaps.mat'],'file')
    disp('Loading pre-calculated permutations ...')
    
    try
        % check if we ran the full (non-reduced) commonality permutations
        load([resultsdir '\permaps.mat']);
    catch
        
        edges_cat = 13:432;
        edges_block = 253:432;
        load([resultsdir '\permaps.mat']);
        
        
        % we couldn't submit all data points because of size limitations, so we
        % are keep only the relevant parts and set all others to zero
        CMEGMRIblock_perm = double(CMEGMRIblock_perm);
        CMEGMRIcat_perm = double(CMEGMRIcat_perm);
        sz = size(CMEGMRIblock_perm); sz(1) = 612;
        % fill in zeros
        CMEGMRIblock_perm = cat(1,zeros(edges_block(1)-1,sz(2),sz(3)), double(CMEGMRIblock_perm));
        CMEGMRIblock_perm = cat(1,CMEGMRIblock_perm, zeros(sz(1)-edges_block(end),sz(2),sz(3)));
        CMEGMRIcat_perm = cat(1,zeros(edges_cat(1)-1,sz(2),sz(3)), double(CMEGMRIcat_perm));
        CMEGMRIcat_perm = cat(1,CMEGMRIcat_perm, zeros(sz(1)-edges_cat(end),sz(2),sz(3)));
        CMEGMRIcue_perm = cue(1,zeros(edges_cue(1)-1,sz(2),sz(3)), double(CMEGMRIcat_perm));
        CMEGMRIcue_perm = cue(1,CMEGMRIcat_perm, zeros(sz(1)-edges_cat(end),sz(2),sz(3)));
        
    end
    
    disp('done.')
    
else
    

    disp('Making 5000 permutations of MEG matrices...')
    RDV_EEG_perm = zeros([size(RDV_EEG) nperms]);
    ind = tril(true(size(RDM_EEG(:,:,1))),-1);
    for i_perm = 1:nperms
        rp = randperm(8);
        curr_RDM = RDM_EEG(rp,rp,:);
        % convert to vector (faster than repeatedly calling mysquareform)
        for i_time = 1:142 % important: use same permutation across time!
            tmp = curr_RDM(:,:,i_time);
            RDV_EEG_perm(:,i_time,i_perm) = tmp(ind);
        end
    end
    
    % get ranks by sorting (because all values are unique in MEG) (this is
    % quite memory intense, alternatively replace the later call to
    % allranks with a call to
    % for i_time = 612:-1:1, for i_perm = 5000:-1:1, allranks(:,i_time,i_perm) = tiedrank2(RDV_MEG_perm(:,i_time,i_perm)); end, end
    disp('calculating ranks in parallel, if you get out of memory, change code in line 90')
    [~,sortind] = sort(RDV_EEG_perm);
    clear RDV_MEG_perm
    [n_cells,n_time] = size(RDV_EEG);
    sortind2 = sortind(:) + kron((0:n_cells:(nperms*n_time*n_cells-1))',ones(n_cells,1));
    clear sortind
    clear allranks
    allranks(sortind2) = repmat((1:n_cells)',n_time*nperms,1);
    clear sortind2
    allranks = reshape(allranks,[n_cells n_time nperms]);
    
    disp('Running 5000 permutations...')
    disp('(come back in a few hours)')
    
    ct = 0;
    ct2 = 0;
    fprintf(repmat('\f',1,9))
    for j_roi = length(roi):-1:1
        
        xMRI = tiedrank2(RDV_fMRI(:,j_roi));
        xcat = tiedrank2(RDV_category);
        xcue = tiedrank2(RDV_cue);
        
        for i_perm = nperms:-1:1
            ct = ct+1;
            if ~mod(ct,5) % five models
                ct2 = ct2+1;
                fprintf(repmat('\b',1,9))
                fprintf('%04i/%04i',ct2,nperms)
            end
            
            for i_time = size(RDM_EEG,3):-1:1
                if i_time >4 && i_time <size(RDM_EEG,3)-4
                y = mean(allranks(:,i_time-4:i_time+4, i_perm),2);
                else
                    y = allranks(:,i_time,i_perm);
                end
               
                
                % and now again with another calculation

                rEEG_MRIcuecat= correlate([y xMRI xcue xcat],'type','spearman','method','semipartialcorr');

                rEEG_MRIcat = correlate([y xMRI xcat],'type','spearman','method','semipartialcorr');
                rEEG_MRIcue = correlate([y xMRI(:,1) xcue], 'type','spearman','method','semipartialcorr');


                CMEGMRIcat_perm(i_time,j_roi, i_perm) = rEEG_MRIcuecat(2,1).^2-rEEG_MRIcue(2,1).^2;
                CMEGMRIcue_perm(i_time,j_roi, i_perm) = rEEG_MRIcuecat(2,1).^2-rEEG_MRIcat(2,1).^2;
%                 
                
                
            end
        end
        
    end
    
    fprintf('\ndone.\n')

end


%% Cluster statistics

% get cluster correction using 0.05 as cluster-inducing threshold and
% restricting our analysis to time periods where we have a hypothesis (see
% below)

for i_roi = length(roi):-1:1
    tmp = squeeze(CMEGMRIcat_perm(:,i_roi,:));
    tmp_sorted = sort(tmp','descend')'; %#ok<*UDIM>
    cutoff_cat(:,i_roi) = tmp_sorted(:,floor(alpha*size(tmp,2)));
    mean_cat(:,i_roi) = tmp_sorted(:,floor(0.5*size(tmp,2)));
    
    
    tmp = squeeze(CMEGMRIcue_perm(:,i_roi,:));
    tmp_sorted = sort(tmp','descend')';
    cutoff_cue(:,i_roi) = tmp_sorted(:,floor(alpha*size(tmp,2)));
    mean_cue(:,i_roi) = tmp_sorted(:,floor(0.5*size(tmp,2)));
    
    % the cutoff below is a relatively conservative estimate (no permutation is larger)
    % cutoff_task(:,i_roi) = tmp_sorted(:,1);
end

% loop over permutations to get maximum cluster size
% pick edges: for task: 13:432 (101ms to 3600ms which is 0 to 3500)
%             for block: 253:432 (2101ms to 3600ms)
edges_cat = 1:142;
edges_block = 1:142;
edges_cue = 1:142;
for i_perm = nperms:-1:1
    for j_roi = length(roi):-1:1
    % find clusters gives us the cluster sizes of the current permutation
    % in the current ROI that are larger than the cutoff value
    c_cat(i_perm,j_roi) = max(find_clusters(CMEGMRIcat_perm(edges_cat,j_roi,i_perm)>cutoff_cat(edges_cat,j_roi)));
    c_cue(i_perm,j_roi) = max(find_clusters(CMEGMRIcue_perm(edges_cue,j_roi,i_perm)>cutoff_cue(edges_cue,j_roi)));
    end
end

for j_roi = length(roi):-1:1
    c_sorted = sort(c_cat(:,j_roi),'descend');
    clust_cutoff_cat(j_roi) = c_sorted(floor(alpha*size(c_cat,1)));
    
    
    c_sorted = sort(c_cue(:,j_roi),'descend');
    clust_cutoff_cue(j_roi) = c_sorted(floor(alpha*size(c_cue,1)));
    
end

% now get cutoff as maximum cutoff of all (which we need for a cutoff corrected for multiple comparisons)
clust_cutoff_catall = max(clust_cutoff_cat);
clust_cutoff_cueall = max(clust_cutoff_cue);

% with these cutoffs check out real cluster sizes
% TODO: we may potentially get percentiles for cluster sizes found empirically for p-values?

all_clustind_cat = zeros(142,length(roi));
all_clustind_cue = zeros(142,length(roi));
for j_roi = length(roi):-1:1
    [c,~,~,clustind] = find_clusters(CMEGMRIcat(edges_cat,j_roi)>cutoff_cat(edges_cat,j_roi));
    for i_c = 1:length(c) % loop over cluster sizes
        if c(i_c)<=clust_cutoff_catall % was <clust_cutoff_cat(j_roi)
            clustind(clustind==i_c) = 0;
        end
    end
    all_clustind_cat(edges_cat,j_roi) = clustind;
    
    
    [c,~,~,clustind] = find_clusters(CMEGMRIcue(edges_cue,j_roi)>cutoff_cue(edges_cue,j_roi));
    for i_c = 1:length(c)
        if c(i_c)<=clust_cutoff_cueall % was <clust_cutoff_block(j_roi)
            clustind(clustind==i_c) = 0;
        end
    end
    all_clustind_cue(edges_cue,j_roi) = clustind;
end

%% Save Results

if save_perms
    save ([resultsdir '\permaps'], 'CMEGMRIcat_perm', 'CMEGMRIcue_perm');     
end    
save ([resultsdir '\results'], 'CMEGMRIcat', 'CMEGMRIcue', 'Fusionmat', 'roi'); 
save ([resultsdir '\stats'], 'all_clustind_cue', 'all_clustind_cat'); 


%% Plotting

% Let's plot results (cat and block within one, i.e. 5 plots)
% using a quadratic scale y-axis and plotting significance above
% (quadratic scale y-axis to reflect an axis similar to correlation values
% of previous RSA and fusion results)

tlim = [-100 1550];
tind = 1:142;

smooth_kern = 3; % smoothed results may look nicer (no smoothing, i.e. kernel = 1 used for figure in paper)
% y = normpdf(linspace(-2.355,2.355,smooth_kern));
y = normpdf(linspace(-2.355,0,smooth_kern));
% normalize to 1
y = y/sum(y);
if length(y)==1
    y = 1;
end

% below is not commented well, but it's only plotting
for i_roi = 1:length(roi)

%     if i_roi == 1 || i_roi == 2 || i_roi == 3
        tytick = [-0.3:0.1:0.8]; % was -0.01:0.01:0.08
        tytick(4) = 0;
%     elseif i_roi == 4
%         tytick = [-0.01 0 0.01 0.04:0.02:0.1 0.15 0.20]; % was 0.04:0.02:0.20
%     elseif i_roi == 5
%         tytick = [-0.01 0 0.01 0.05:0.05:0.2 0.3 0.4]; % was 0.05:0.05:0.35
%     end
%     
    tyticklabel = num2cell(tytick);
    tylim = tytick([1 end]);
    
    CMEGMRIcat (t<0,:) = CMEGMRIcat (t<0,:)./2;
    catplot = CMEGMRIcat(tind,i_roi);
    
    CMEGMRIcue (t<0,:) = CMEGMRIcue (t<0,:)./2;
    cueplot = CMEGMRIcue(tind,i_roi);
    
    Fusionmat (:,t<0) = Fusionmat (:,t<0)./2;
    ceilplot = abs(Fusionmat(i_roi,tind).^2);%.^2; %i_roi, tind
%     ceilplot = rEEG_MRI(tind, i_roi)';
    
    catplot = [catplot((smooth_kern-1)/2:-1:1); catplot; catplot(end:-1:end-(smooth_kern-1)/2+1)]; %#ok<AGROW>
    catplot = conv(catplot,y,'valid');
    
    
    
    cueplot = [cueplot((smooth_kern-1)/2:-1:1); cueplot; cueplot(end:-1:end-(smooth_kern-1)/2+1)]; %#ok<AGROW>
    cueplot = conv(cueplot,y,'valid');
    
    ceilplot = [ceilplot((smooth_kern-1)/2:-1:1), ceilplot, ceilplot(end:-1:end-(smooth_kern-1)/2+1)]; %#ok<AGROW>
    ceilplot = conv(ceilplot,y,'valid');
    ceilplot = ceilplot';
    
    hf = figure;
    h3 = area(t(1:142),sqrt(ceilplot)); % plot fusion as baseline

    hold on
    % we convert values to quadratic scale and change the y-axis
    s = sign(catplot);
    h1 = plot(t(1:142),s.*sqrt(abs(catplot)));

    
    
    s = sign(cueplot);
    h2 = plot(t(1:142),s.*sqrt(abs(cueplot)));
    xlim(tlim)
%     legend({'fusion','cat','blockect'})
    title(strrep(roi{i_roi},'_','\_'))
    h1.LineWidth = 2.5;
    h1.Color = [177,38,32]/255; %red 

    h2.LineWidth = 2.5;
    h2.Color = [230,210,0]/255; %yellow
    h3.LineWidth = 1.5;
    h3.FaceColor = [0.9 0.9 0.9];
    h3.EdgeColor = [0 0 0];
    ha = gca;
    ha.YTickMode = 'manual';
    ha.YTick = tytick;
    %ha.YTick = sign(tytick).*sqrt(abs(tytick));
    ha.YTickLabel = tyticklabel;
    ha.YLim = tylim;
   % ha.YLim = sign(tylim).*sqrt(abs(tylim));
    
    % add significance bars
    sig_cat = double(all_clustind_cat(tind,i_roi)>0);
    sig_cat(sig_cat==0) = NaN; % NaN will allow plotting as a line with breaks

    sig_cue = double(all_clustind_cue(tind,i_roi)>0);
    sig_cue(sig_cue==0) = NaN;

    % fact will provide where to place the significance bars
    fact_cat = [-0.07 -0.07 -0.07 -0.07 -0.07 -0.07]; ; % was 0.075 0.075 0.075 0.19 0.335
    fact_cue  = [-0.05 -0.05 -0.05 -0.05 -0.05 -0.05]; % was 0.07 0.07 0.07 0.18 0.32
    
    ht1 = plot(t,sig_cat*fact_cat(i_roi),'color', [177,38,32]/255,'linewidth',1.5);
    ht2 = plot(t,sig_cue*fact_cue(i_roi),'color', [230,210,0]/255,'linewidth',1.5);
    
    % plot time 0
    plot([0 0],ha.YLim,'k--','linewidth',1.5)
    plot(xlim,[0 0],'k--','linewidth',1.5)
    set(gcf,'color', 'w')
    set(gcf,'Position', [300,300, 750, 500])
    
% activate below if you want to save the results    
%     print(['commonality_' roinames{i_roi} '.eps'],'-painters','-depsc')
% saveas(gcf,['commonality_' roinames{i_roi} '.png'],'png')
    

end
