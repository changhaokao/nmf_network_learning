function merge_subject_corrMat(sublist, window_setting)

% sublist
if nargin<1 | isempty(sublist)
    sublist = textread('sublist', '%s');
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end
nSub = numel(sublist);

% window setting
if nargin<2
    window_width = 10; % unit: TR
    window_overlap = 8;
else
    window_width = window_setting.window_width;
    window_overlap = window_setting.window_overlap;
end


% analysis type
analysisType = 'pearson';

% directory
dirType = 'data_24mc_CSF_WM';
dirCorr = fullfile(dirType, 'corr_pearson');
dirGroup = fullfile(dirType, 'corr_pearson_group');
mkdir(dirGroup);


% setting
nRun = 4;
nTR_used = 220;
nWindow = (nTR_used-window_width)/(window_width-window_overlap) + 1;
nWindow_total = nWindow*nRun;

ROI_unused = [2, 4, 5, 8, 10, 78, 81, 83, 120, 179, 180, 181, 182, 247, 248, 249, 250];
nROI = 264-numel(ROI_unused);
nCorr = nROI*(nROI-1)/2;

corrMat_group = zeros(nCorr,nWindow_total*nSub);


for s = 1:nSub
    
    % subname
    subname = sublist{s};
    
    % load file
    filename = fullfile(dirCorr, sprintf('%s_corr_%s_window_%d_overlap_%d.mat', subname, analysisType, window_width, window_overlap));
    load(filename);
    
    % merge data
    col_range = [(s-1)*nWindow_total+1:s*nWindow_total];
    corrMat_group(:,col_range) = corrMat;
    
end

% output
outputFile = fullfile(dirGroup, sprintf('group_corr_%s_window_%d_overlap_%d.mat', analysisType, window_width, window_overlap));
save(outputFile, 'corrMat_group', '-v7.3');

% positive
corrMat_group_pos = corrMat_group;
corrMat_group_pos((corrMat_group<0)) = 0;
outputFile = fullfile(dirGroup, sprintf('group_pos_corr_%s_window_%d_overlap_%d.mat', analysisType, window_width, window_overlap));
save(outputFile, 'corrMat_group_pos', '-v7.3');

% negative
corrMat_group_neg = corrMat_group;
corrMat_group_neg((corrMat_group>0)) = 0;
outputFile = fullfile(dirGroup, sprintf('group_neg_corr_%s_window_%d_overlap_%d.mat', analysisType, window_width, window_overlap));
save(outputFile, 'corrMat_group_neg', '-v7.3');

clear corrMag_group;
