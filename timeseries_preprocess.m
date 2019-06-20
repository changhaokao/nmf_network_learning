function timeseries_preprocess(sublist)

% regress out the influence of motion parameters and signals from CSF and WM

% sublist
if nargin<1 | isempty(sublist)
    sublist = textread('sublist','%s');
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end
nSub = numel(sublist);


% directory
dirTS = fullfile('timeseries', 'diameter_9mm');
dirSeg = fullfile('timeseries', 'CSF_WM_mask_erode');
dirMC = 'mc';
dirReg = 'reg';
mkdir(dirReg);


% setting
nRun = 4;

samplerate = 1/2.5;
freq_BOLD = [0.01, 0.08];
freq_BOLD = (freq_BOLD/samplerate)*2;
n_order = 1;


% preprocess data
for s = 1:nSub
    tic

    % subname
    subname = sublist{s};
    fprintf('%s\n',subname);
    
    % load file
    filename = sprintf('%s_ts.mat',subname);
    load(fullfile(dirTS, filename));
    
    % load confound files
    mc_all = cell(1,nRun);
    CSF_all = cell(1,nRun);
    WM_all = cell(1,nRun);
    for r = 1:nRun
        
        %%%%% load data %%%%%
        mc_file = sprintf('%s_r%.2d_st_mc.par', subname, r);
        mc_original = load(fullfile(dirMC, mc_file));
        mc_original_sq = mc_original.^2;
        mc_pre = [zeros(1,6);mc_original(1:end-1,:)];
        mc_pre_sq = mc_pre.^2;
        mc_par = [mc_original, mc_original_sq, mc_pre, mc_pre_sq]; % mc_{t}, mc_{t}^2, mc_{t-1}, mc_{t-1}^2

        seg_file = sprintf('%s_ts_seg.mat',subname);
        load(fullfile(dirSeg, seg_file)); % ts_seg.CSF, ts_seg.WM
        
        %%%%% filter %%%%%
        [b,a] = butter(n_order, freq_BOLD, 'bandpass');
        
        mc_par_filtered = filtfilt(b,a,mc_par);
        mc_all{r} = mc_par_filtered;

        CSF_filtered = filtfilt(b,a,ts_seg.CSF.ts_mean{r});
        CSF_all{r} = CSF_filtered;

        WM_filtered = filtfilt(b,a,ts_seg.WM.ts_mean{r});
        WM_all{r} = WM_filtered;
        
    end

    
    % ROI
    nROI = numel(ts.ROI);
    clear ts_reg
    for i = 1:nROI
        for r = 1:nRun
            
            % data
            ts_raw = ts.ROI(i).ts_mean{r};
            
            if ~any(isnan(ts_raw))
                
                % filter
                [b,a] = butter(n_order, freq_BOLD, 'bandpass');
                ts_filtered = filtfilt(b,a,ts_raw);
                
                % regress out motion parameter
                x = [mc_all{r}, CSF_all{r}, WM_all{r}];
                y = ts_filtered;
                beta = glmfit(x,y,'normal');
                
                yhat = glmval(beta,x,'identity');
                y_reg = y-yhat;
                
                % save result
                ts_reg.ROI(i).ts_mean{1,r} = y_reg;
                ts_reg.ROI(i).isNaN(1,r) = 0;
                
            else
                
                ts_reg.ROI(i).ts_mean{1,r} = nan(size(ts_raw));
                ts_reg.ROI(i).isNaN(1,r) = 1;
                
            end
            
        end
    end
    
    
    % save variable
    outputfile = sprintf('%s_ts_reg', subname);
    save(fullfile(dirReg,outputfile),'ts_reg');
    
    toc
end

