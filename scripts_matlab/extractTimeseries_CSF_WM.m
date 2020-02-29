function extractTimeseries_CSF_WM(sublist)

% directory
dataDir = fullfile('preprocessed');
segDir = fullfile('structure_segmentation');

tsDir = fullfile('timeseries',sprintf('CSF_WM_mask_erode'));
mkdir(tsDir);

% subject list
if nargin<1
    sublist = textread('sublist', '%s');    
else
    if ~iscell(sublist)
        sublist = {sublist};
    end
end
nSub = numel(sublist);

% numer of run
nRun = 4;

for s = 1:nSub
    tic
    % subject
    subname = sublist{s};
    subDir = fullfile(dataDir,subname);

    % show progress
    fprintf('%s\n', subname);
    

    %%%%% load CSF & WM %%%%%
    % CSF
    filename = sprintf('%s_CSF_mask_erode.nii', subname);
    imgFile = fullfile(segDir, subname, filename);
    V = spm_vol(imgFile);
    Voxel2MNI = V(1).mat;
    mat = spm_read_vols(V);
    nZ = size(mat,3);
    xall = [];
    yall = [];
    zall = [];
    for z = 1:nZ
        [x,y] = find(mat(:,:,z)==1);
        xall = [xall;x];
        yall = [yall;y];
        zall = [zall;ones(size(x))*z];
    end
    voxel_coordinate = [xall';yall';zall'];
    MNI_coordinate = Voxel2MNI*[voxel_coordinate; ones(1,size(voxel_coordinate,2))];
    MNI_CSF = MNI_coordinate(1:3,:);

    % WM
    filename = sprintf('%s_WM_mask_erode.nii', subname);
    imgFile = fullfile(segDir, subname, filename);
    V = spm_vol(imgFile);
    Voxel2MNI = V(1).mat;
    mat = spm_read_vols(V);
    nZ = size(mat,3);
    xall = [];
    yall = [];
    zall = [];
    for z = 1:nZ
        [x,y] = find(mat(:,:,z)==1);
        xall = [xall;x];
        yall = [yall;y];
        zall = [zall;ones(size(x))*z];
    end
    voxel_coordinate = [xall';yall';zall'];
    MNI_coordinate = Voxel2MNI*[voxel_coordinate; ones(1,size(voxel_coordinate,2))];
    MNI_WM = MNI_coordinate(1:3,:);
    
    %%%%% extract timeseries %%%%%
    for r = 1:nRun

        % run
        runname = sprintf('r%.2d', r);
        
        % load image file
        filename = sprintf('%s_%s_st_mc_nLin_sm6_ds_in.nii', subname, runname);
        imgFile = fullfile(subDir,filename);
        
        % get MNI2Voxel matrix
        V = spm_vol(imgFile);
        MNI2Voxel = inv(V(1).mat);
        
          
        % timeseries: CSF
        MNI_extract = MNI_CSF;
        voxel_extract = MNI2Voxel*[MNI_extract; ones(1, size(MNI_extract,2))];
        voxel_extract = round(voxel_extract);
        voxel_extract = voxel_extract(1:3,:);
        val = spm_get_data(imgFile, voxel_extract);
        
        ts_seg.CSF.ts_mean{r} = average_ts(val);

        % timeseries: WM
        MNI_extract = MNI_WM;
        voxel_extract = MNI2Voxel*[MNI_extract; ones(1, size(MNI_extract,2))];
        voxel_extract = round(voxel_extract);
        voxel_extract = voxel_extract(1:3,:);
        val = spm_get_data(imgFile, voxel_extract);
        
        ts_seg.WM.ts_mean{r} = average_ts(val);
            
    end
    
    tsFile = fullfile(tsDir, sprintf('%s_ts_seg', subname));
    save(tsFile, 'ts_seg');
    toc
end


function ts_mean = average_ts(ts)

nCol = size(ts,2);
idx_include = ones(nCol,1);

for i = 1:nCol
    
    if all(ts(:,i)==0) | any(isnan(ts(:,i)))
        
        idx_include(i) = 0;
        
    end
    
end
idx_include = logical(idx_include);
ts = ts(:,idx_include);
ts_mean = mean(ts,2);




