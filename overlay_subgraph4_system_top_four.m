function overlay_subgraph4_system_top_four()

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirSubsystem = fullfile(dirType, 'result_subsystem_10_8');

% setting
dim = [61, 73, 61]; % 3mm X 3mm X 3mm

% method of ROI
ROI_radius = 100; % mm

% directory
% ROI list
ROI_filename = 'ROI.mat';
load(ROI_filename);

% list_system
list_system = [6,8,9,13];
nSystem = numel(list_system);

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
                    brain_node(i,j,k) = ROI(idx_ROI).subsystem_idx;
                end
                
            end
            
        end
    end
end


matrix_val = zeros(dim);
for i = 1:nSystem
    
    system_number = list_system(i);
    idx_system = (brain_node==system_number);
    matrix_val(idx_system) = i;
    
end


% output file
outputfile = fullfile(dirSubsystem, 'node_subgraph_4_systemn_top_four.nii');
V.fname = outputfile;
V.private.dat.fname = V.fname;
spm_write_vol(V, matrix_val);

gzip(outputfile);
delete(outputfile);




