function merge_pos_neg(window_setting)

dirType = 'data_24mc_CSF_WM';

if nargin<1
    window_width = 10; % unit: TR
    window_overlap = 8;
else
    window_width = window_setting.window_width;
    window_overlap = window_setting.window_overlap;
end

dirGroup = fullfile(dirType, 'corr_pearson_group');

filename = sprintf('group_pos_corr_pearson_window_%d_overlap_%d.mat', window_width, window_overlap);
load(fullfile(dirGroup, filename));

filename = sprintf('group_neg_corr_pearson_window_%d_overlap_%d.mat', window_width, window_overlap); 
load(fullfile(dirGroup, filename));
    
corrMat_group_both = abs([corrMat_group_pos,corrMat_group_neg]);

outputFile = sprintf('group_both_corr_pearson_window_%d_overlap_%d.mat', window_width, window_overlap);
save(fullfile(dirGroup, outputFile), 'corrMat_group_both', '-v7.3');
