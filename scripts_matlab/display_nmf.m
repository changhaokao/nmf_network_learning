function display_nmf()

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirFig = fullfile(dirType, 'figures_display_10_8');
dirConsensus = fullfile(dirType, 'result_consensus_10_8');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% interaction X time %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = fullfile(dirType, 'group_both_corr_pearson_window_10_overlap_8.mat');
load(filename);

[nEdge, nTime] = size(corrMat_group_both);

% graph
figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*2;
fg.position(4) = fg.position(4)*1.5;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
imagesc(corrMat_group_both(:,[1:100:27000]),[0,1]);
colormap('parula');
xlabel('Time windows', 'fontsize', 40);
ylabel('Edges', 'fontsize', 40);
set(gca,...
    'YDir','normal',...
    'XTick',[],...
    'YTick',[],...
    'linewidth',2,...
    'box','off');

% save figure
outputfile = fullfile(dirFig, sprintf('nmf_edgeXtime'));
print(outputfile,'-depsc','-opengl');


nROI = 247;
for idx_example = 1:6
    matrix_subnet = NaN(nROI,nROI);
    current_corr_edge = corrMat_group_both(:,100*idx_example);
    idx = 0;
    for i = 1:nROI-1
        for j = i+1:nROI
            
            idx = idx + 1;
            matrix_subnet(i,j) = current_corr_edge(idx);
            
        end
    end
    matrix_subnet = matrix_subnet';
    figure;
    fg = fig_setting_default();
    set(gca,'YDir','reverse');
    hold on
    
    imagesc(matrix_subnet,[0,1]);
    colormap('parula');
    
    fill([0.5;nROI+0.5;nROI+0.5],[0.5;0.5;nROI+0.5],'w');

    plot([0.5;nROI+0.5;nROI+0.5],[0.5;0.5;nROI+0.5],...
        'color', 'w',...
        'linewidth', 5);
    plot([0.5;0.5;nROI+0.5],[0.5;nROI+0.5;nROI+0.5],...
        'color', [0.7,0.7,0.7],...
        'linewidth', 5);
    plot([0;nROI+1],[0;nROI+1],...
        'color', [0.7,0.7,0.7],...
        'linewidth', 5);
    
    
    hold off
    
    set(gca,...
        'XTick',[],...
        'YTick',[],...
        'linewidth',2,...
        'box','off');
    axis off
    axis image
    xlim([0.5,nROI+0.5]);
    ylim([0.5,nROI+0.5]);
    camroll(-45)
    
    % save figure
    outputfile = fullfile(dirFig, sprintf('example_corr_map_%d', idx_example));
    print(outputfile,'-depsc','-opengl');
end




figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*0.5;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
axis off
colormap('parula');
c = colorbar('fontsize',40, 'Ticks', [0:0.2:1]);
c.Label.String = 'Edge strength';
c.LineWidth = 2;
c.Position(1) = 0.5;
c.Label.Position = [-3,0.5,0];

% save figure
outputfile = fullfile(dirFig, 'nmf_colorbar_corr');
print(outputfile,'-depsc','-opengl');





%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% edge X subgraph %%%%%
%%%%% subgraph X time %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_subgraph = 4;
filename = fullfile(dirConsensus,'result_ens.mat');
load(filename);

nSubgraph = size(subnet_ens,1);
nEdge = size(subnet_ens,2);

% edge X subgraph
% graph
figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*0.5;
fg.position(4) = fg.position(4)*1.5;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);

hold on

imagesc(subnet_ens',[0,0.02]);
xlabel('Subgraphs', 'fontsize', 40);
ylabel('Edges', 'fontsize', 40);
colormap('parula');

plot([-0.5;-0.5;0.5;0.5]+idx_subgraph,[0;nEdge;nEdge;0],...
    'color', 'k',...
    'linewidth', 5);

hold off
set(gca,...
    'YDir','normal',...
    'XTick',[],...
    'YTick',[],...
    'linewidth',2,...
    'box','on');
xlim([0.5,nSubgraph+0.5]);
ylim([0,nEdge]);

% save figure
outputfile = fullfile(dirFig, sprintf('nmf_edgeXsubgraph'));

print(outputfile,'-depsc','-opengl');



% subgraph X time
% graph
figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*1.25;
fg.position(4) = fg.position(4)*0.85;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
hold on
for i = 1:nSubgraph
    
    current_data = coef_ens(i,[1:500]);
    current_data = rescale_data(current_data);
    current_data = current_data + i;
    plot(current_data,...
        'color', [0.4,0.4,0.4],...
        'linewidth', 2);
    
end
hold off

ylim([0,nSubgraph+1]);
xlabel('Time windows', 'fontsize', 40);
ylabel('Subgraphs', 'fontsize', 40);

set(gca,...
    'YDir','normal',...
    'XTick',[],...
    'YTick',[],...
    'linewidth',2);

% save figure
outputfile = fullfile(dirFig, sprintf('nmf_subgraphXtime'));
print(outputfile,'-depsc','-opengl');




%%%%%%%%%%%%%%%%%%%%%%
%%%%% timeseries %%%%%
%%%%%%%%%%%%%%%%%%%%%%
filename = fullfile(dirType, 'bb914_ts_reg.mat');
load(filename);

nROI = 247;
nROI_display = 20;

figure;
fig_setting_default();
hold on
for i = 1:nROI_display
    
    current_data = [ts_reg.ROI(i).ts_mean{1};ts_reg.ROI(i).ts_mean{2};ts_reg.ROI(i).ts_mean{3};ts_reg.ROI(i).ts_mean{4}];
    current_data = current_data(1:200,:);
    current_data = rescale_data(current_data);
    
    current_data = current_data + i;
    plot(current_data,...
        'color', [0.4,0.4,0.4],...
        'linewidth', 2);
    
end
hold off

ylim([0,nROI_display+1]);
xlabel('Time', 'fontsize', 32);
ylabel('BOLD magnitude', 'fontsize', 32);
set(gca,...
    'YDir','normal',...
    'XTick',[],...
    'YTick',[],...
    'linewidth',2);

% save figure
outputfile = fullfile(dirFig, sprintf('nmf_timeseries'));
print(outputfile,'-depsc','-opengl');





function data = rescale_data(data)

data = (data-min(data))./(max(data)-min(data));


