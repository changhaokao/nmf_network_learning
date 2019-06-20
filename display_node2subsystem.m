function display_node2subsystem()

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirSubsystem = fullfile(dirType, 'result_subsystem_10_8');
dirFig = fullfile(dirType, 'figures_display_10_8');
mkdir(dirFig);

cmap_system = [
     0     0     0
     0    73    73
   255   182   119
    73     0   146
     0   109   219
   182   109   255
   109   182   255
   182   219   255
   146     0     0
   146    73     0
   219   209     0
    36   255    36
   255   255   109
    ];
cmap_system = cmap_system./255;
list_system = {
    'Uncertain';
    'Sensory';
    'Cingulo-opercular';
    'Auditory';
    'Default mode';
    'Memory retrieval';
    'Visual';
    'Fronto-parietal';
    'Salience';
    'Subcortical';
    'Ventral attention';
    'Cerebellar';
    'Dorsal attention'
    };
nSystem = numel(list_system);

s = 4;
%%%%% node %%%%%
% load file
filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node.mat',s));
load(filename);

val_min = min(matrix_subnet(:));
val_max = max(matrix_subnet(:));
matrix_subnet = (matrix_subnet-val_min)./(val_max-val_min);

figure;
fig_setting_default;

imagesc(matrix_subnet,[0,1]);
colormap('parula');

set(gca,...
    'XTick',[60:60:240],...
    'YTick',[60:60:240],...
    'linewidth',2,...
    'box','on',...
    'fontsize', 20);
xlabel('Regions', 'fontsize', 20);
ylabel('Regions', 'fontsize', 20);

axis image
pos_node = get(gca,'position');

outputfile = fullfile(dirFig, sprintf('connectivity_node_subgraph%.2d',s));
print(outputfile,'-depsc','-opengl');


% colorbar
figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*0.3;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
imagesc(matrix_subsystem_idx(:,1,1));
colormap(cmap_system);
pos_system = get(gca,'Position');
pos_system([2,4]) = pos_node([2,4]);
pos_system([3]) = 0.1;
set(gca,'Position',pos_system);

list_tick = NaN(nSystem,1);
for i = 1:nSystem
    idx = find(matrix_subsystem_idx(:,1,1)==i);
    list_tick(i) = mean(idx);
end

for i = 1:nSystem
    text(1.8,list_tick(i),list_system{i},'fontsize',12);
end

axis off

outputfile = fullfile(dirFig, 'connectivity_subsystem_index');
print(outputfile,'-depsc','-opengl');







