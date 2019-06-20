function extractTimeseries()

% directory
dataDir = fullfile('preprocessed');
tsDir = fullfile('timeseries','diameter_9mm');
mkdir(tsDir);

% subject list
sublist = textread('sublist', '%s');
nSub = numel(sublist);

% ROI list
ROI_filename = 'ROI.mat';
load(ROI_filename);
nROI = numel(ROI);

% numer of run
nRun = 4;

% method of ROI
setting.shape = 'sphere';
setting.diameter = 3; % 3 voxles (3mm X 3mm X 3mm), i.e., 9mm based on Power et al., 2011, Neuron.
ROI_template = create_template(setting);
nCoordinate = size(ROI_template, 2);

for s = 1:nSub
    
    % subject
    subname = sublist{s};
    subDir = fullfile(dataDir,subname);
    
    for r = 1:nRun
        
        % run
        runname = sprintf('r%.2d', r);
        
        % load image file
        filename = sprintf('%s_%s_st_mc_nLin_sm6_ds_in.nii', subname, runname);
        imgFile = fullfile(subDir,filename);
        
        % get MNI2Voxel matrix
        V = spm_vol(imgFile);
        MNI2Voxel = inv(V(1).mat);
        
        for i = 1:nROI
            
            % show progress
            fprintf('%s, %s, %d\n', subname, runname, i);
            
            % create coordinates to extract timeseries
            MNI_center = reshape(ROI(i).MNI_coordiante, 3, 1);
            voxel_center = MNI2Voxel*[MNI_center; ones(1, 1)];
            voxel_center = voxel_center(1:3,:);
            voxel_extract = ROI_template + repmat(voxel_center, 1, nCoordinate);
            voxel_extract = round(voxel_extract);
            
            % extract timeseries
            val = spm_get_data(imgFile, voxel_extract);
            
            % save to mat
            ts.ROI(i).MNI_center = MNI_center;
            ts.ROI(i).voxel_center = voxel_center;
            ts.ROI(i).voxel_extract = voxel_extract;
            ts.ROI(i).ts_each{r} = val;
            ts.ROI(i).ts_mean{r} = average_ts(val);
            
        end
    end
    
    tsFile = fullfile(tsDir, sprintf('%s_ts', subname));
    save(tsFile, 'ts');
    
end




function ROI_template = create_template(setting)

nNeighbor = (setting.diameter-1)/2;
vector = [-nNeighbor:nNeighbor];
n = numel(vector);
cube_voxels = [sort(repmat(vector, 1, n^2)); repmat(sort(repmat(vector, 1, n)), 1, n); repmat(vector, 1, n^2)];

switch lower(setting.shape)
    case 'sphere'
        distance = (sum(cube_voxels.^2, 1)).^(1/2);
        idx = (distance<=nNeighbor);
        x = cube_voxels(:, idx);
    case 'cube'
        x = cube_voxels;
end
ROI_template = x;

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




