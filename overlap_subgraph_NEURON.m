function allData = overlap_subgraph_NEURON()

window_width = 10;
window_overlap = 8;

nPerm = 10000;

% directory
data_type = 'result_raw';
dirType = fullfile('result_24mc_CSF_WM', data_type);
dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));
dirFig = fullfile(dirType, sprintf('figures_%d_%d', window_width, window_overlap));
mkdir(dirFig);

% list system-connectoin
list_system = {'all'};
nSystem = numel(list_system);

% list_NEURON
list_NEURON = {
    'Fig3A_CPP_unthresh';
    'Fig3B_RU_unthresh';
    'Fig3C_Rwd_unthresh';
    };
nType = numel(list_NEURON);

% ROI list
ROI_filename = 'ROI.mat';
load(ROI_filename);

% sample
samplefile = fullfile('NEURON-D-14-01017_data','Fig3A_CPP_thresh.nii.gz');
sampleimg = load_nii(samplefile);
voxel2MNI = [sampleimg.hdr.hist.srow_x;sampleimg.hdr.hist.srow_y;sampleimg.hdr.hist.srow_z;[0,0,0,1]];
MNI2voxel = inv(voxel2MNI);

% merge MNI coordinate
idx_include = [ROI.idx_include];
nValid = sum(idx_include);
ROI = ROI(idx_include==1);
MNI_all = NaN(nValid, 3);
for i = 1:nValid
    if ROI(i).idx_include==1
        MNI_all(i,:) = ROI(i).MNI_coordinate;
    end
end
voxel_all = MNI2voxel*[MNI_all,ones(nValid,1)]';
voxel_all = voxel_all(1:3,:)';
voxel_all = round(voxel_all);


% subgraph
s = 4;
% s = 8;


% load file
filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_node',s));
load(filename);

% within, between
idx_between = (matrix_subsystem_idx(:,:,1)~=matrix_subsystem_idx(:,:,2));
matrix_subnet_within = matrix_subnet;
matrix_subnet_within(idx_between) = nan;

idx_within = (matrix_subsystem_idx(:,:,1)==matrix_subsystem_idx(:,:,2));
matrix_subnet_between = matrix_subnet;
matrix_subnet_between(idx_within) = nan;

% all
node_strength.all = nansum(matrix_subnet,1)/2;
val_min = min(node_strength.all);
val_max = max(node_strength.all);
node_strength.all = (node_strength.all-val_min)./(val_max-val_min);

% within
node_strength.within = nansum(matrix_subnet_within,1)/2;
val_min = min(node_strength.within);
val_max = max(node_strength.within);
node_strength.within = (node_strength.within-val_min)./(val_max-val_min);

% between
node_strength.between = nansum(matrix_subnet_between,1)/2;
val_min = min(node_strength.between);
val_max = max(node_strength.between);
node_strength.between = (node_strength.between-val_min)./(val_max-val_min);

nROI = numel(ROI_idx);

for t = 1:nType
    
    typename = list_NEURON{t};
    filename = fullfile('NEURON-D-14-01017_data',sprintf('%s.nii.gz', typename));
    cope = load_nii(filename);
    data = cope.img;
    data = mean(data,4);
    
    val_fmri = NaN(nROI,1);
    val_subgraph.all = NaN(nROI,1);
    val_subgraph.within = NaN(nROI,1);
    val_subgraph.between = NaN(nROI,1);
    for r = 1:nROI
        
        idx_roi = (ROI_idx==ROI(r).idx);
        current_voxel = voxel_all(r,:);
        
        x = current_voxel(1);
        y = current_voxel(2);
        z = current_voxel(3);
        val_fmri(r) = data(x,y,z);
        val_subgraph.all(r) = node_strength.all(idx_roi);
        val_subgraph.within(r) = node_strength.within(idx_roi);
        val_subgraph.between(r) = node_strength.between(idx_roi);
        
        
        
    end
    
    allData.(typename).activation = val_fmri;
    allData.(typename).connectivity = val_subgraph;
    
end




% graph
for c = 1:nSystem
    
    systemname = list_system{c};
    
    for t = 1:nType
        
        typename = list_NEURON{t};
        
        figure;
        fg = fig_setting_default();
        
        x = allData.(typename).connectivity.(systemname);
        y = allData.(typename).activation;
        [rho, p] = permutation_corr(x,y,nPerm);
        xval = linspace(min(x),max(x),100)';
        [param, yhat, ci] = polypredci(x, y, 1, 0.95, xval);
        yhat = double(yhat);
        ci = double(ci);
        
        hold on
        [hl, hp] = boundedline(xval, yhat, ci, 'alpha');
        
        set(hl,...
            'Linesmoothing', 'on',...
            'Color', 'r',...
            'LineStyle', '-',...
            'LineWidth', 4);
        
        set(hp,...
            'Linesmoothing', 'on',...
            'FaceColor', 'r');
        
        plot(x,y,...
            'linestyle', 'none',...
            'linewidth', 1,...
            'Marker', 'o',...
            'MarkerSize', 8,...
            'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', [0.5,0.5,0.5]);
        
        hold off
        
        set(gca, 'fontsize', 28.57, 'linewidth', 2);
        
        xlabel('Normalized edge strength', 'fontsize', 28.57);
        ylabel('Activation (z-statistic)', 'fontsize', 28.57);
        
        switch typename
            case {'Fig3A_CPP_unthresh'}
                title(sprintf('CPP'), 'fontsize',28.57);
                xrange = [-0.1,1.1];
                yrange = ylim;
                xlim(xrange);
                ylim(yrange);
                xpos = xrange(1)+(xrange(2)-xrange(1))*0.75;
                ypos = yrange(1)+(yrange(2)-yrange(1))*0.15;
            case {'Fig3B_RU_unthresh'}
                title(sprintf('RU'), 'fontsize',28.57);
                xrange = [-0.1,1.1];
                yrange = ylim;
                xlim(xrange);
                ylim(yrange);
                xpos = xrange(1)+(xrange(2)-xrange(1))*0.75;
                ypos = yrange(1)+(yrange(2)-yrange(1))*0.15;
            case {'Fig3C_Rwd_unthresh'}
                title(sprintf('Reward'), 'fontsize',28.57);
                xrange = [-0.1,1.1];
                %             yrange = [-0.1,1.1];
                yrange = ylim;
                xlim(xrange);
                ylim(yrange);
                xpos = xrange(1)+(xrange(2)-xrange(1))*0.05;
                ypos = yrange(1)+(yrange(2)-yrange(1))*0.9;
        end
        
        set(gca, 'XTick', [0:0.2:1]);
        
        if p>=.001
            text_p = sprintf('p=%.3f',p);
        elseif p<.001 & p>=.0001
            text_p = 'p<.001';
        elseif p<.0001
            text_p = 'p<.0001';
        end
        
        text(xpos,ypos,sprintf('r=%.3f\n%s', rho, text_p), 'fontsize',28.57);
        
        outputfile = fullfile(dirFig, sprintf('overlap_NEURON_%s_subgraph%d_%s', systemname, s, typename));
        print(outputfile,'-depsc','-opengl');
        
    end
end




