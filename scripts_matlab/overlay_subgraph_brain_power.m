function overlay_subgraph_brain_power()

% directory
dirSubsystem = 'result_subsystem_Power';

% setting
dim = [61, 73, 61]; % 3mm X 3mm X 3mm

% method of ROI
ROI_radius = 4.5; % mm

% directory
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
brain_ROI = zeros(dim);
brain_system = zeros(dim);
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
                    brain_ROI(i,j,k) = ROI(idx_ROI).idx;
                    brain_system(i,j,k) = ROI(idx_ROI).subsystem_idx;
                end
                
            end
            
        end
    end
end

% output file
outputfile = fullfile(dirSubsystem, 'power_brain_ROI.nii');
V.fname = outputfile;
V.private.dat.fname = V.fname;
spm_write_vol(V, brain_ROI);
gzip(outputfile);
delete(outputfile);

outputfile = fullfile(dirSubsystem, 'power_brain_system.nii');
V.fname = outputfile;
V.private.dat.fname = V.fname;
spm_write_vol(V, brain_system);
gzip(outputfile);
delete(outputfile);





