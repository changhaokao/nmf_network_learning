function run_subsystem_permutation(window_setting)

if nargin<1
    window_width = 10;
    window_overlap = 8;
else
    window_width = window_setting.window_width;
    window_overlap = window_setting.window_overlap;
end

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
nSubnet = 10;

% dirType = 'result_24mc_CSF_WM/result_glm_bold_temporal_derivative';
% nSubnet = 9;

dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));

% setting
nSubsystem = 13;
nPermutation = 10000;

% sort subsystem
tic
for s = 1:nSubnet
    
    fprintf('subgraph %.2d:\n',s);
    
    %%%%% load data %%%%%
    % subgraph: node
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node.mat', s));
    load(filename)
    nNode = size(matrix_subnet,1);
    
    % subgraph: subsystem
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_subsystem.mat', s));
    load(filename)
    
    subsystem_interaction = subsystem_subnet;
    subsystem_strength = sum(subsystem_subnet);
    
    
    %%%%% permutation %%%%%
    subsystem_interaction_perm = NaN(nSubsystem,nSubsystem,nPermutation);
    subsystem_strength_perm = NaN(nPermutation, nSubsystem);
    
    for p = 1:nPermutation
       
        % permute node index
        idx_perm = randperm(nNode);
        current_subsystem_idx = matrix_subsystem_idx;
        current_subsystem_idx(:,:,1) = current_subsystem_idx(idx_perm,:,1);
        current_subsystem_idx(:,:,2) = current_subsystem_idx(:,idx_perm,2);
        
        % assign nodes to subsystems
        subsystem_perm = zeros(nSubsystem,nSubsystem);
        for i = 1:nSubsystem
            for j = 1:nSubsystem
                
                idx = (current_subsystem_idx(:,:,1)==i)&(current_subsystem_idx(:,:,2)==j);
                current_data = matrix_subnet(idx);
                subsystem_perm(i,j) = mean(current_data);
                
            end
        end
        
        subsystem_interaction_perm(:,:,p) = subsystem_perm;
        subsystem_strength_perm(p,:) = sum(subsystem_perm,1);
        
    end
    
    %%%%% permutation test: subsystem %%%%%
    perm_ttest = (subsystem_strength_perm>repmat(subsystem_strength,nPermutation,1));
    pval = (mean(perm_ttest));
    fprintf('sys%.2d\t',[1:nSubsystem]);
    fprintf('\n');
    fprintf('%.4f\t',pval);
    fprintf('\n');
    
    % save data
    outputfile = fullfile(dirSubsystem, sprintf('subsystem_perm_subgraph_%.2d.mat',s));
    save(outputfile, 'subsystem_strength_perm', 'subsystem_strength', 'pval');
    
    
    %%%%% permutation test: interaction %%%%%
    perm_ttest = (subsystem_interaction_perm>repmat(subsystem_interaction,[1,1,nPermutation]));
    pval = (mean(perm_ttest,3));
    
    % save data
    outputfile = fullfile(dirSubsystem, sprintf('subsystem_interaction_perm_subgraph_%.2d.mat',s));
    save(outputfile, 'subsystem_interaction_perm', 'subsystem_interaction', 'pval');
    
end
toc



