function display_Power()

cfgfile = 'option_Power.mat';
load(cfgfile);
modular = EC.nod.ModularNumber;
n = numel(modular);

color_subsystem = load('color_subsystem');
color = color_subsystem(modular,[2:4]);
color = color/255;
EC.nod.CMm([1:n]',:) = color;
EC.nod.CMm_temp([1:n]',:) = color;
EC.nod.CM([1:n]',:) = color;

EC.nod.size = 1;
EC.nod.size_size = 3;

EC.msh.alpha = 0.4;

save(cfgfile,'EC');

nodefile = 'node_Power.node';
BrainNet_MapCfg('BrainMesh_ICBM152.nv',nodefile,cfgfile);

outputfile = 'fig_Power_node';
% export_fig(outputfile,'-dpng');
print(outputfile,'-depsc','-opengl');




