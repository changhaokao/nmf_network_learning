window_setting.window_width = 10;
window_setting.window_overlap = 8;

result_type = 'result_raw';
% result_type = 'result_glm_bold_temporal_derivative';
idx_fig = 1;
typename = '';

%%%%%%%%%%%%%%%%%%%%%
%%%%% subsystem %%%%%
%%%%%%%%%%%%%%%%%%%%%
identify_subsystem_connectivity(window_setting); % identify subsystem from nodes
run_subsystem_permutation(window_setting); % run permutation to identify significant system-by-system connection
run_subsystem_boostrapping(window_setting); % run boostraping to estimate confidence interval of system-by-system connection
group_subsystem_connectivity(window_setting); % summarize subgraph as node-by-node and system-by-system connection

create_subgraph_node_list(); % create node list for significant system-by-system connection
display_brain_connectivity(); % display nodes for significant system-by-system connection by BrainNet



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% task-related effect %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r-square between within-subject and between-subject effects
group_rsquared_within_between_effect();

% the relationship between trial-by-trial task variables and subgraph
% expression and the relationship between individual normative learning and
% subgraph expression
group_behav_subgraph_coef(typename, idx_fig, result_type, window_setting);

% contributions of within-system edges and between-system edges to
% trial-by-trial task effects and individual normative learning; also,
% contributions of each functional system were investigated as well.
group_behav_subgraph_coef_comparison();

% contributions of eash system-by-systems to task effects and individual
% normative learning
group_behav_subgraph_coef_remove_interaction();



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% map a subgraph on the brain %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overlay_subgraph_brain(); % create brain image of edge strength
run_surfice_cmd_NEURON(); % create brain map from McGuire et al. (2014; Neuron) by surfice
run_surfice_cmd_subgraph(); % create brain map of the target subgraph by surfice
run_clip_fmri_figures(); % merge figures created by surfice
overlap_subgraph_NEURON(); % relationship between edge strength and activation effects


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% method visualization %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_nmf(); % NMF
display_node2subsystem(); % node-by-node connection to system-by-system connection


