function overlay_subgraph_brain()

window_width = 10;
window_overlap = 8;

% setting
dim = [61, 73, 61]; % 3mm X 3mm X 3mm

nSubgraph = 10;

% method of ROI
ROI_radius = 100; % mm

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
% dirType = 'result_24mc_CSF_WM/result_residual';
dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));

% ROI list
ROI_filename = 'ROI.mat';
load(ROI_filename);

% mask
maskfile = 'MNI152_T1_3mm_brain_mask.nii.gz';
V = spm_vol(maskfile);
matrix_mask = spm_read_vols(V);

% sample
samplefile = fullfile('NEURON-D-14-01017_data','Fig3A_CPP_thresh.nii.gz');
V = spm_vol(samplefile);
voxel2MNI = V.mat;
MNI2voxel = inv(voxel2MNI);

% merge MNI coordinate
idx_include = [ROI.idx_include];
nValid = sum(idx_include);
ROI = ROI(idx_include==1);
MNI_all = NaN(nValid, 3);
for i = 1:nValid
    if ROI(i).idx_include==1
        MNI_all(i,:) = ROI(i).MNI_coordinate;
    end
end

% assign node index based on distance
brain_node = zeros(dim);
for i = 1:dim(1)
    for j = 1:dim(2)
        for k = 1:dim(3)
            
            if matrix_mask(i,j,k)==1
                
                current_voxel = [i;j;k;1];
                current_MNI = voxel2MNI*current_voxel;
                current_MNI = current_MNI(1:3)';
                current_MNI = repmat(current_MNI,nValid,1);
                distance_node = sqrt(sum((MNI_all-current_MNI).^2,2));
                [val,idx] = min(distance_node);
                
                if val<=ROI_radius
                    idx_ROI = idx(1);
                    brain_node(i,j,k) = idx_ROI;
                end
                
            end
            
        end
    end
end

% subgraph: pos
for s = 1:nSubgraph
    
    matrix_val = zeros(dim);
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node',s));
    load(filename);
    
    node_strength = sum(matrix_subnet,1)/2;
    val_min = min(node_strength);
    val_max = max(node_strength);
    node_strength = (node_strength-val_min)./(val_max-val_min);
    
    nROI = numel(ROI_idx);
    
    for r = 1:nROI
        
        idx_roi = (ROI_idx==ROI(r).idx);
        idx_voxel = (brain_node==r);
        
        matrix_val(idx_voxel) = node_strength(idx_roi);
        
    end
    
    outputfile = fullfile(dirSubsystem,sprintf('node_strength_subgraph%.2d_pos.nii',s));
    V.fname = outputfile;
    V.private.dat.fname = V.fname;
    spm_write_vol(V, matrix_val);
    
    gzip(outputfile);
    delete(outputfile);
    
end





% subgraph: neg
for s = 1:nSubgraph
    
    matrix_val = zeros(dim);
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node',s));
    load(filename);
    
    node_strength = sum(matrix_subnet,1)/2;
    val_min = min(node_strength);
    val_max = max(node_strength);
    node_strength = (node_strength-val_min)./(val_max-val_min);
    node_strength = abs(node_strength-1);
    
    nROI = numel(ROI_idx);
    
    for r = 1:nROI
        
        idx_roi = (ROI_idx==ROI(r).idx);
        idx_voxel = (brain_node==r);
        
        matrix_val(idx_voxel) = node_strength(idx_roi);
        
    end
    
    outputfile = fullfile(dirSubsystem,sprintf('node_strength_subgraph%.2d_neg.nii',s));
    V.fname = outputfile;
    V.private.dat.fname = V.fname;
    spm_write_vol(V, matrix_val);
    
    gzip(outputfile);
    delete(outputfile);
    
end

