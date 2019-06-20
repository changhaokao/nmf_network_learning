function corrTimeseries(sublist, window_setting)

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

tsDir = fullfile(dirType, 'reg');

dirCorr = fullfile(dirType, 'corr_pearson');
mkdir(dirCorr);


% parameter
nRun = 4;
ROI_unused = [2, 4, 5, 8, 10, 78, 81, 83, 120, 179, 180, 181, 182, 247, 248, 249, 250];
TR_unused = [1:6];
nTR_unused = numel(TR_unused);
TR = 2.5; % second

% Band-pass filter for wavelet coherence analysis
if strcmp(analysisType, 'wavelet')
    bandPass = [0.02, 0.08];
end


for s = 1:nSub
    
    tic
    % subject
    subname = sublist{s};
    fprintf('%d, %s\n',s, subname);
    filename = sprintf('%s_ts_reg.mat',subname);
    load(fullfile(tsDir, filename));
    ts = ts_reg;
    
    % ROI
    ts.ROI(ROI_unused) = [];
    nROI = numel(ts.ROI);
    nTR_total = numel(ts.ROI(1).ts_mean{1});
    nTR_used = nTR_total-nTR_unused;
    
    % initialize corr matrix
    nCorr = nROI*(nROI-1)/2;
    nWindow = (nTR_used-window_width)/(window_width-window_overlap) + 1;
    
    % merge timeseries from all ROI
    for r = 1:nRun
        mergeROI{r} = zeros(nTR_used, nROI);
        for idx_roi = 1:nROI
            tsROI = ts.ROI(idx_roi).ts_mean{r};
            tsROI(TR_unused) = [];
            mergeROI{r}(:,idx_roi) = tsROI;
        end
    end
    
    
    corrMat = zeros(nCorr, nWindow*nRun);
    
    for r = 1:nRun
        for w = 1:nWindow
            
            idx_time = (r-1)*nWindow + w;
            
            % period to calculate corr
            periodTR = [1:window_width] + (window_width-window_overlap)*(w-1);
            
            % Pearson's correlation
            periodData = mergeROI{r}(periodTR, :);
            rho = corr(periodData);
            
            idx_row = 0;
            for i = 1:nROI-1
                for j = i+1:nROI
                    
                    idx_row = idx_row + 1;
                    corrMat(idx_row,idx_time) = rho(i,j);
                    
                end
            end
            
        end % end of window
    end % end of run
    
    outputFile = fullfile(dirCorr, sprintf('%s_corr_%s_window_%d_overlap_%d.mat', subname, analysisType, window_width, window_overlap));
    save(outputFile, 'corrMat');
    
    toc
    
end % end of subject


