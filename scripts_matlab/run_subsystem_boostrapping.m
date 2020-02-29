function run_subsystem_boostrapping(window_setting)

if nargin<1
    window_width = 10;
    window_overlap = 8;
else
    window_width = window_setting.window_width;
    window_overlap = window_setting.window_overlap;
end

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
nSubgraph = 10;

% dirType = 'result_24mc_CSF_WM/result_glm_bold_temporal_derivative';
% nSubgraph = 9;

dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));

% setting
nBoostrap = 10000;

list_system = {
    '1. Uncertain';
    '2. Sensory';
    '3. Cingulo-opercular';
    '4. Auditory';
    '5. Default mode';
    '6. Memory retrieval';
    '7. Visual';
    '8. Fronto-parietal';
    '9. Salience';
    '10. Subcortical';
    '11. Ventral attention';
    '12. Cerebellar';
    '13. Dorsal attention';
    };
nSystem = numel(list_system);



% start
within_connection_boostrap = NaN(nSubgraph,nBoostrap);
between_connection_boostrap = NaN(nSubgraph,nBoostrap);

system_connection_boostrap = NaN(nSubgraph,nSystem,nBoostrap);
system_within_connection_boostrap = NaN(nSubgraph,nSystem,nBoostrap);
system_between_connection_boostrap = NaN(nSubgraph,nSystem,nBoostrap);
tic
for s = 1:nSubgraph
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node.mat',s));
    load(filename);
    
    %%%%% within-connection %%%%%
    idx_edge = (matrix_subsystem_idx(:,:,1)==matrix_subsystem_idx(:,:,2));
    edge_val = matrix_subnet(idx_edge);
    nEdge = numel(edge_val);
    for i = 1:nBoostrap
        
        idx_random = ceil(rand(nEdge,1)*nEdge);
        val_mean = mean(edge_val(idx_random));
        within_connection_boostrap(s,i) = val_mean;
        
    end
    
    %%%%% between-connection %%%%%
    idx_edge = (matrix_subsystem_idx(:,:,1)~=matrix_subsystem_idx(:,:,2));
    edge_val = matrix_subnet(idx_edge);
    nEdge = numel(edge_val);
    for i = 1:nBoostrap
        
        idx_random = ceil(rand(nEdge,1)*nEdge);
        val_mean = mean(edge_val(idx_random));
        between_connection_boostrap(s,i) = val_mean;
        
    end
    
    %%%%% system %%%%%
    for k = 1:nSystem
        
        % all
        idx_edge = (matrix_subsystem_idx(:,:,1)==k);
        edge_val = matrix_subnet(idx_edge);
        nEdge = numel(edge_val);
        for i = 1:nBoostrap
            
            idx_random = ceil(rand(nEdge,1)*nEdge);
            val_mean = mean(edge_val(idx_random));
            system_connection_boostrap(s,k,i) = val_mean;
            
        end
        
        % within-connection
        idx_edge = (matrix_subsystem_idx(:,:,1)==k & matrix_subsystem_idx(:,:,2)==k);
        edge_val = matrix_subnet(idx_edge);
        nEdge = numel(edge_val);
        for i = 1:nBoostrap
            
            idx_random = ceil(rand(nEdge,1)*nEdge);
            val_mean = mean(edge_val(idx_random));
            system_within_connection_boostrap(s,k,i) = val_mean;
            
        end
        
        % between-connection
        idx_edge = (matrix_subsystem_idx(:,:,1)==k & matrix_subsystem_idx(:,:,2)~=k);
        edge_val = matrix_subnet(idx_edge);
        nEdge = numel(edge_val);
        for i = 1:nBoostrap
            
            idx_random = ceil(rand(nEdge,1)*nEdge);
            val_mean = mean(edge_val(idx_random));
            system_between_connection_boostrap(s,k,i) = val_mean;
            
        end
        
    end
    
end

toc


% output file
outputfile = fullfile(dirSubsystem, 'subsystem_boostrap.mat');
save(outputfile,...
    'within_connection_boostrap', 'between_connection_boostrap',...
    'system_connection_boostrap', 'system_within_connection_boostrap', 'system_between_connection_boostrap');





