function modify_options(option_setting)

view_direction = option_setting.view_direction;
nSubgraph = option_setting.nSubgraph;
dirSubsystem = option_setting.dirSubsystem;

% if nargin<1
%     p_threshold = .99;
% end

opacity = 0.4;

color_subsystem = load('result_subsystem_Power/color_subsystem');

% subsystem: interaction
for i = 1:nSubgraph
    
    cfgfile = fullfile(dirSubsystem, sprintf('option_node_interaction_subgraph_%.2d.mat',i));
    load(cfgfile);
    modular = EC.nod.ModularNumber;
    n = numel(modular);
    
    color = color_subsystem(modular,[2:4]);
    color = color/255;
    EC.nod.color = 3;
    EC.nod.CMm([1:n]',:) = color;
    EC.nod.CMm_temp([1:n]',:) = color;
    EC.nod.CM([1:n]',:) = color;
    
    EC.nod.size = 1;
    EC.nod.size_size = 3;
    
    
    EC.edg.draw = 2;
    
    EC.lot.view_direction = view_direction;
    
    EC.msh.alpha = opacity;
    
    save(cfgfile,'EC');
    
end


