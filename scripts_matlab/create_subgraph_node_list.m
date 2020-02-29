function create_subgraph_node_list(window_setting)

if nargin<1
    window_width = 10;
    window_overlap = 8;
else
    window_width = window_setting.window_width;
    window_overlap = window_setting.window_overlap;
end

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));
nSubgraph = 10;

% setting
size_node = 3;

% ROI
load('ROI.mat');
ROI_unused = [2, 4, 5, 8, 10, 78, 81, 83, 120, 179, 180, 181, 182, 247, 248, 249, 250];
ROI(ROI_unused) = [];
nROI = numel(ROI);

for s = 1:nSubgraph
    
    %%%%% subsystem %%%%%
    p_threshold = 0.05/13;
    % load file
    filename = fullfile(dirSubsystem, sprintf('subsystem_perm_subgraph_%.2d',s));
    load(filename);
    
    % open output file
    outputfile = fullfile(dirSubsystem, sprintf('node_subgraph_%.2d.node',s));
    fid = fopen(outputfile, 'w');
    
    % assign nodes
    for r = 1:nROI
        
        subsystem_idx = ROI(r).subsystem_idx;
        MNI_coordinate = ROI(r).MNI_coordinate;
        
        idx_ROI = ROI(r).idx_include;
        idx_pval = (pval(subsystem_idx)<p_threshold);
        idx_include = (idx_ROI&idx_pval);
        
        
        if idx_include
            fprintf(fid, '%.2f\t', MNI_coordinate);
            fprintf(fid, '%d\t', subsystem_idx);
            fprintf(fid, '%d\t', size_node);
            fprintf(fid, '-\n');
        end
        
    end
    fclose(fid);
    
    
    %%%%% subsystem: interaction %%%%%
    p_threshold = 0.05/91;
    % load file
    filename = fullfile(dirSubsystem, sprintf('subsystem_interaction_perm_subgraph_%.2d',s));
    load(filename);
    
    % load interaction
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node.mat',s));
    load(filename)
    
    % open output file: node
    outputfile = fullfile(dirSubsystem, sprintf('node_interaction_subgraph_%.2d.node',s));
    fid = fopen(outputfile, 'w');
    
    subsystem_sig = zeros(size(pval,1),1);
    [x,y] = find(pval<p_threshold);
    subsystem_sig(unique([x;y])) = 1;
    subsystem_sig = logical(subsystem_sig);
    
    % assign nodes
    idx_node = [];
    for r = 1:nROI
        
        subsystem_idx = ROI(r).subsystem_idx;
        MNI_coordinate = ROI(r).MNI_coordinate;
        
        idx_ROI = ROI(r).idx_include;
        idx_pval = (subsystem_sig(subsystem_idx)==1);
        idx_include = (idx_ROI&idx_pval);
        
        if idx_include
            idx_node = [idx_node, find(ROI_idx==r)];
            fprintf(fid, '%.2f\t', MNI_coordinate);
            fprintf(fid, '%d\t', subsystem_idx);
            fprintf(fid, '%d\t', size_node);
            fprintf(fid, '-\n');
        end
        
    end
    fclose(fid);
    
    
    % assign edges
    matrix_edge = matrix_subnet(idx_node,idx_node);
    
    % open output file: edges
    outputfile = fullfile(dirSubsystem, sprintf('edge_interaction_subgraph_%.2d.edge',s));
    dlmwrite(outputfile, matrix_edge,'delimiter','\t');
    
    
    
end


