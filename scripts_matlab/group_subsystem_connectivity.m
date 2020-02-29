function group_subsystem_connectivity(window_setting)

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
dirFig = fullfile(dirType, sprintf('figures_%d_%d', window_width, window_overlap));

mkdir(dirFig);


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

% setting
clims = [0,0.01];

%%%%% node %%%%%
for s = 1:nSubgraph
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node.mat',s));
    load(filename);
    
    %     min(min(matrix_subnet))
    %     max(max(matrix_subnet))
    
    figure;
    fig_setting_default;
    %     subplot(3,3,idx);
    %     imagesc(matrix_subnet,clims);
    imagesc(matrix_subnet,[0,0.02]);
    colormap('parula');
    %     colorbar;
    
    %     title(sprintf('subgraph%.2d: node',s));
    set(gca,...
        'XTick',[],...
        'YTick',[],...
        'linewidth',2,...
        'box','on');
    
    axis image
    
    outputfile = fullfile(dirFig, sprintf('connectivity_node_subgraph%.2d',s));
    %     export_fig(outputfile, '-png');
    print(outputfile,'-depsc','-opengl');
    
end

%%%% subsystem %%%%%
for s = 1:nSubgraph
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_subsystem.mat',s));
    load(filename);
    val_min = min(subsystem_subnet(:));
    val_max = max(subsystem_subnet(:));
    subsystem_subnet = (subsystem_subnet-val_min)./(val_max-val_min);
    
    figure;
    fig_setting_default;
    imagesc(subsystem_subnet, [0,1]);
    colormap('parula');
    
    set(gca,...
        'XTick',[1:2:13],...
        'YTick',[1:2:13],...
        'fontsize',20,...
        'linewidth',2,...
        'box','on');
    xlabel('Systems', 'fontsize', 20);
    ylabel('Systems', 'fontsize', 20);
    
    axis image
    
    %     title(sprintf('subgraph%.2d: subsystem',s));
    
    outputfile = fullfile(dirFig, sprintf('connectivity_subsystem_subgraph%.2d',s));
    print(outputfile,'-depsc', '-opengl');
    
end


%%%%% subsystem: threshold by interaction %%%%%
p_threshold = 0.05/91;
for s = 1:nSubgraph
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_subsystem.mat',s));
    load(filename);
    val_min = min(subsystem_subnet(:));
    val_max = max(subsystem_subnet(:));
    subsystem_subnet = (subsystem_subnet-val_min)./(val_max-val_min);
    
    filename = fullfile(dirSubsystem, sprintf('subsystem_interaction_perm_subgraph_%.2d.mat',s));
    load(filename);
    idx_threshold = (pval<p_threshold);
    subsystem_subnet(~idx_threshold) = NaN;
    
    figure;
    fg = fig_setting_default;
    
    imagesc(subsystem_subnet, [0,1]);
    custom_color = parula;
    custom_color(1,:) = [1,1,1];
    colormap(custom_color);
    
    set(gca,...
        'XTick',[1:2:13],...
        'YTick',[1:2:13],...
        'fontsize',30,...
        'linewidth',2,...
        'box','on');
    
    title(sprintf('Subgraph %.1d',s),'fontsize',33.33);
    xlabel('Systems', 'fontsize', 33.33);
    ylabel('Systems', 'fontsize', 33.33);
    
    axis image
    
    
    outputfile = fullfile(dirFig, sprintf('connectivity_subsystem_threshold_subgraph%.2d',s));
    print(outputfile,'-depsc','-opengl');
    
end


% figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*0.5;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
axis off
custom_color = parula;
custom_color(1,:) = [1,1,1];
colormap(custom_color);
c = colorbar('fontsize',30, 'Ticks', [0:0.2:1]);
c.Label.String = 'Normalized edge strength';
c.LineWidth = 2;
c.Position(1) = 0.5;
c.Label.Position = [-2,0.5,0];
% save figure
outputfile = fullfile(dirFig, 'connectivity_subsystem_colorbar_threshold');
print(outputfile,'-depsc','-opengl');


% figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*0.5;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
axis off
custom_color = parula;
colormap(custom_color);
c = colorbar('fontsize',30, 'Ticks', [0:0.2:1]);
c.Label.String = 'Normalized edge strength';
c.LineWidth = 2;
c.Position(1) = 0.5;
c.Label.Position = [-2,0.5,0];
% save figure
outputfile = fullfile(dirFig, 'connectivity_subsystem_colorbar');
print(outputfile,'-depsc','-opengl');






%%%%% subsystem: within-connection & between-connection %%%%%
within_connection = NaN(nSubgraph,1);
between_connection = NaN(nSubgraph,1);

system_connection = NaN(nSubgraph,nSystem);
system_within_connection = NaN(nSubgraph,nSystem);
system_between_connection = NaN(nSubgraph,nSystem);
for s = 1:nSubgraph
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_subsystem.mat',s));
    load(filename);
    
    system_connection(s,:) = mean(subsystem_subnet,1);
    
    val_diag = diag(subsystem_subnet);
    within_connection(s) = mean(val_diag);
    system_within_connection(s,:) = val_diag;
    
    val_off_diag = triu(subsystem_subnet,1);
    idx_include = logical(triu(ones(size(subsystem_subnet)),1));
    val_off_diag = val_off_diag(idx_include);
    between_connection(s) = mean(val_off_diag);
    
    system_between_connection(s,:) = (sum(subsystem_subnet,1) - val_diag')/12;
    
    
end
within_ratio = (within_connection-between_connection)./(within_connection+between_connection);


%%%%% boostrap %%%%%
filename = fullfile(dirSubsystem, sprintf('subsystem_boostrap'));
load(filename);

within_ratio_boostrap = (within_connection_boostrap-between_connection_boostrap)./(within_connection_boostrap+between_connection_boostrap);



%%%%% graph %%%%%
facecolor = [0.5,0.5,0.5];
edgecolor = 'none';
errorbarcolor = 'k';
barwidth = 0.7;


% within_ratio
barData = within_ratio;
[barData, idx_subgraph] = sort(barData,'ascend');

semData = std(within_ratio_boostrap,0,2);
semData = semData(idx_subgraph,1);

figure;
fig_setting_default;

hold on

h = barh(barData);
set(h,...
    'facecolor', facecolor,...
    'edgecolor', edgecolor,...
    'linewidth', 2,...
    'barwidth', barwidth);

barsx = h.XData;
e = errorbar(barData, barsx, semData, 'horizontal');
set(e,...
    'color', errorbarcolor,...
    'linestyle', 'none',...
    'linewidth', 2);

plot([0;0],[0;nSubgraph+0.5],...
    'color','k',...
    'linestyle','-',...
    'linewidth',2);

hold off


for j = 1:nSubgraph
    if barData(j)>0
        xpos = -0.06;
    elseif barData(j)<0
        xpos = 0.06;
    end
    text(xpos, j, num2str(idx_subgraph(j)),...
        'HorizontalAlignment', 'Center',...
        'color', 'k',...
        'fontsize', 24,...
        'background','w');
end
xlim([-0.7,0.7]);
ylim([0.5, nSubgraph+0.5]);
set(gca,'YColor','w');
set(gca,'XTick',[-0.6,-0.4,-0.2,0,0.2,0.4,0.6]);
set(gca,'linewidth',2, 'fontsize', 24);
xlabel('Relative strength','fontsize',24);

outputfile = fullfile(dirFig, 'relative_strength');
print(outputfile,'-depsc','-opengl');



% within_connection
barData = within_connection;
[barData, idx_subgraph] = sort(barData,'ascend');

semData = std(within_connection_boostrap,0,2);
semData = semData(idx_subgraph,1);

figure;
fig_setting_default;

hold on
h = barh(barData);
set(h,...
    'facecolor', facecolor,...
    'edgecolor', edgecolor,...
    'linewidth', 2,...
    'barwidth', barwidth);

barsx = h.XData;
e = errorbar(barData, barsx, semData, 'horizontal');
set(e,...
    'color', errorbarcolor,...
    'linestyle', 'none',...
    'linewidth', 2);
hold off

xlim([0,0.01]);
ylim([0.5, nSubgraph+0.5]);
set(gca,'XTick', [0:0.005:0.01], 'XTickLabel', [0:0.005:0.01]);
set(gca,'YTick', [1:nSubgraph], 'YTickLabel', idx_subgraph);
set(gca,'linewidth',2, 'fontsize', 24);
xlabel('Within-system strength','fontsize',24);

outputfile = fullfile(dirFig, 'within_system_strength');
print(outputfile,'-depsc','-opengl');



% between_connection
barData = between_connection;
[barData, idx_subgraph] = sort(barData,'ascend');

semData = std(between_connection_boostrap,0,2);
semData = semData(idx_subgraph,1);


figure;
fig_setting_default;

hold on
h = barh(barData);
set(h,...
    'facecolor', facecolor,...
    'edgecolor', edgecolor,...
    'linewidth', 2,...
    'barwidth', barwidth);

barsx = h.XData;
e = errorbar(barData, barsx, semData, 'horizontal');
set(e,...
    'color', errorbarcolor,...
    'linestyle', 'none',...
    'linewidth', 2);
hold off

xlim([0,0.006]);
ylim([0.5, nSubgraph+0.5]);
set(gca,'XTick', [0:0.002:0.006], 'XTickLabel', [0:0.002:0.006]);
set(gca,'YTick', [1:nSubgraph], 'YTickLabel', idx_subgraph);
set(gca,'linewidth',2, 'fontsize', 24);
xlabel('Between-system strength','fontsize',24);

outputfile = fullfile(dirFig, 'between_system_strength');
print(outputfile,'-depsc','-opengl');




%%%% system %%%%%
idx_subgraph = 4;

% system_connection
barData = system_connection(idx_subgraph,:);
[barData, idx_system] = sort(barData,'ascend');

current_data = system_connection_boostrap(idx_subgraph,:,:);
current_data = reshape(current_data, nSystem, []);
semData = std(current_data,0,2);
semData = semData(idx_system,1);


figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(4) = fg.position(4)*1.2;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);

hold on
h = barh(barData);
set(h,...
    'facecolor', facecolor,...
    'edgecolor', edgecolor,...
    'linewidth', 2,...
    'barwidth', barwidth);

barsx = h.XData;
e = errorbar(barData, barsx, semData, 'horizontal');
set(e,...
    'color', errorbarcolor,...
    'linestyle', 'none',...
    'linewidth', 2);
hold off


xlim([0,0.01]);
ylim([0.5, nSystem+0.5]);
set(gca,'XTick', [0:0.004:0.008], 'XTickLabel', {'0', '0.004', '0.008'});
set(gca,'YTick', [1:nSystem], 'YTickLabel', list_system(idx_system));
set(gca,'linewidth',2, 'fontsize', 22);
xlabel('Average system strength','fontsize',33.33);

outputfile = fullfile(dirFig, sprintf('system_strength_%.2d', idx_subgraph));
print(outputfile,'-depsc');

% system_within_connection
barData = system_within_connection(idx_subgraph,:);
[barData, idx_system] = sort(barData,'ascend');

current_data = system_within_connection_boostrap(idx_subgraph,:,:);
current_data = reshape(current_data, nSystem, []);
semData = std(current_data,0,2);
semData = semData(idx_system,1);

figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(4) = fg.position(4)*1.2;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);

hold on
h = barh(barData);
set(h,...
    'facecolor', facecolor,...
    'edgecolor', edgecolor,...
    'linewidth', 2,...
    'barwidth', barwidth);

barsx = h.XData;
e = errorbar(barData, barsx, semData, 'horizontal');
set(e,...
    'color', errorbarcolor,...
    'linestyle', 'none',...
    'linewidth', 2);
hold off

xlim([0,0.02]);
ylim([0.5, nSystem+0.5]);
set(gca,'XTick', [0:0.01:0.02], 'XTickLabel', [0:0.01:0.02]);
set(gca,'YTick', [1:nSystem], 'YTickLabel', list_system(idx_system));
set(gca,'linewidth',2, 'fontsize', 22);
xlabel('Within-system strength','fontsize',33.33);

outputfile = fullfile(dirFig, sprintf('system_within_system_strength_%.2d', idx_subgraph));
print(outputfile,'-depsc','-opengl');


% system_between_connection
barData = system_between_connection(idx_subgraph,:);
[barData, idx_system] = sort(barData,'ascend');

current_data = system_between_connection_boostrap(idx_subgraph,:,:);
current_data = reshape(current_data, nSystem, []);
semData = std(current_data,0,2);
semData = semData(idx_system,1);

figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(4) = fg.position(4)*1.2;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);

hold on
h = barh(barData);
set(h,...
    'facecolor', facecolor,...
    'edgecolor', edgecolor,...
    'linewidth', 2,...
    'barwidth', barwidth);

barsx = h.XData;
e = errorbar(barData, barsx, semData, 'horizontal');
set(e,...
    'color', errorbarcolor,...
    'linestyle', 'none',...
    'linewidth', 2);
hold off

xlim([0,0.01]);
ylim([0.5, nSystem+0.5]);
set(gca,'XTick', [0:0.004:0.008], 'XTickLabel', {'0', '0.004', '0.008'});
set(gca,'YTick', [1:nSystem], 'YTickLabel', list_system(idx_system));
set(gca,'linewidth',2, 'fontsize', 22);
xlabel('Between-system strength','fontsize',33.33);

outputfile = fullfile(dirFig, sprintf('system_between_system_strength_%.2d', idx_subgraph));
print(outputfile,'-depsc','-opengl');



