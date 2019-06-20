function identify_subsystem_connectivity(window_setting)

if nargin<1
    window_width = 10;
    window_overlap = 8;
else
    window_width = window_setting.window_width;
    window_overlap = window_setting.window_overlap;
end

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
% dirType = 'result_24mc_CSF_WM/result_glm_bold_temporal_derivative';

dirConsensus = fullfile(dirType, sprintf('result_consensus_%d_%d', window_width, window_overlap));
dirResult = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));
mkdir(dirResult);

% load ROI
load('ROI.mat');

% load consensus matrix
load(fullfile(dirConsensus, 'result_ens.mat'));

% setting
ROI_idx = [1:numel(ROI)]';
ROI_unused = [2, 4, 5, 8, 10, 78, 81, 83, 120, 179, 180, 181, 182, 247, 248, 249, 250];
ROI(ROI_unused) = [];
nROI = numel(ROI);
ROI_idx(ROI_unused) = [];

nSubnet = size(subnet_ens,1);


% subsystem
subsystem_idx = [ROI.subsystem_idx]';
nSubsystem = max(subsystem_idx);

[subsystem_idx, original_location] = sort(subsystem_idx);
ROI_idx = ROI_idx(original_location);

% sort subsystem
for s = 1:nSubnet
    
    current_subnet = subnet_ens(s,:);
    
    %%%%% subnet: nodes %%%%%
    % reshape
    matrix_subnet = zeros(nROI,nROI);
    idx = 0;
    for i = 1:nROI-1
        for j = i+1:nROI
            
            idx = idx + 1;
            matrix_subnet(i,j) = current_subnet(idx);
            
        end
    end
    matrix_subnet = matrix_subnet + matrix_subnet';
    
    % sort by subsystem
    original_subnet = matrix_subnet;
    matrix_subnet = zeros(nROI,nROI);
    matrix_subsystem_idx = zeros(nROI,nROI,2);
    for i = 1:nROI-1
        for j = i+1:nROI
            
            row = original_location(i);
            col = original_location(j);
            matrix_subnet(i,j) = original_subnet(row,col);
            
        end
    end
    matrix_subsystem_idx(:,:,1) = repmat(subsystem_idx,1,nROI);
    matrix_subsystem_idx(:,:,2) = matrix_subsystem_idx(:,:,1)';
    
    matrix_subnet = matrix_subnet + matrix_subnet';
    
    % save
    filename = fullfile(dirResult, sprintf('subgraph_%.2d_subnet_node.mat',s));
    save(filename, 'matrix_subnet', 'matrix_subsystem_idx', 'ROI_idx');
    
    
    %%%%% subnet: subsystem %%%%%
    subsystem_subnet = zeros(nSubsystem,nSubsystem);
    for i = 1:nSubsystem
        for j = 1:nSubsystem
            
            idx = (matrix_subsystem_idx(:,:,1)==i)&(matrix_subsystem_idx(:,:,2)==j);
            current_data = matrix_subnet(idx);
            subsystem_subnet(i,j) = mean(current_data);
            
        end
    end
    
    % save
    filename = fullfile(dirResult, sprintf('subgraph_%.2d_subnet_subsystem.mat',s));
    save(filename, 'subsystem_subnet');
    
end




