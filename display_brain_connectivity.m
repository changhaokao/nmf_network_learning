function display_brain_connectivity()

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirSubsystem = fullfile(dirType, 'result_subsystem_10_8');

% list view
list_view = {'sagittal','axial','coronal'};

% setting
p_threshold = 0.99;
view_direction = 1;
nSubgraph = 10;

option_setting.p_threshold = p_threshold;
option_setting.view_direction = view_direction;
option_setting.nSubgraph = nSubgraph;
option_setting.dirSubsystem = dirSubsystem;

modify_options(option_setting);

% subsystem: interaction
for i = 1:nSubgraph
    
    nodefile = fullfile(dirSubsystem, sprintf('node_interaction_subgraph_%.2d.node',i));
    cfgfile = fullfile(dirSubsystem, sprintf('option_node_interaction_subgraph_%.2d.mat',i));
    
    BrainNet_MapCfg('BrainMesh_ICBM152.nv',nodefile,cfgfile);
    
    outputfile = fullfile(dirSubsystem, sprintf('fig_interaction_subgraph_%.2d_%s',i,list_view{view_direction}));
    print(outputfile,'-depsc','-opengl');
    
end

