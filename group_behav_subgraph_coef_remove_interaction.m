function group_behav_subgraph_coef_remove_interaction()


window_setting.window_width = 10;
window_setting.window_overlap = 8;

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_setting.window_width, window_setting.window_overlap));
dirFig = fullfile(dirType, sprintf('figures_%d_%d', window_setting.window_width, window_setting.window_overlap));
mkdir(dirFig);

% setting
nPermutation = 10000;
list_coef = {'rel'};
nCoef = numel(list_coef);

nInteraction = 91;


% load data
fprintf('loading data...\n');
tic
allData.all = group_behav_subgraph_coef('', 0);
for idx_interaction = 1:nInteraction
    allData.lesion{idx_interaction} = group_behav_subgraph_coef('_remove_interaction_individual', 0, 'result_raw', window_setting, idx_interaction);
end
toc

nSub = size(allData.all.reg.coef_behav.cpp_ru_vals.rel,1);
nSystem = 13;

idx_subgraph = 4;

filename = fullfile(dirSubsystem, sprintf('subsystem_interaction_perm_subgraph_%.2d.mat',idx_subgraph));
load(filename);
p_threshold = 0.05/91;
subnet_thresholded = (pval<p_threshold);
idx_edge = logical(triu(ones(size(subnet_thresholded))));
nEdge_threshold = sum(subnet_thresholded(idx_edge));

% color: parula
custom_parula = parula;

% color: red
nLevel = 64;
color_top = [0.5, 0, 0];
color_middle = [1, 0, 0];
color_bottom = [1, 1, 1];
custom_red = NaN(nLevel, 3);
for i = 1:3
    custom_red(1:nLevel/2, i) = linspace(color_top(i), color_middle(i), nLevel/2);
end
for i = 1:3
    custom_red(nLevel/2+1:nLevel, i) = linspace(color_middle(i), color_bottom(i), nLevel/2);
end

% color: blue
nLevel = 64;
color_top = [0, 0, 0.5];
color_middle = [0, 0.5, 1];
color_bottom = [1, 1, 1];
custom_blue = NaN(nLevel, 3);
for i = 1:3
    custom_blue(1:nLevel/2, i) = linspace(color_top(i), color_middle(i), nLevel/2);
end
for i = 1:3
    custom_blue(nLevel/2+1:nLevel, i) = linspace(color_middle(i), color_bottom(i), nLevel/2);
end

% color: redblue
% pval_mc = 0.05/91;
pval_mc = 0.0005;
list_pval = [0.05, 0.01, 0.001, pval_mc];
custom_redblue = [];
nLevel = numel(list_pval);
idx_color = round(linspace(1, 64, nLevel+1));
idx_color = idx_color(1:nLevel);
custom_redblue = [custom_redblue; custom_blue(idx_color, :)];
custom_redblue = [custom_redblue; [1,1,1]];
custom_redblue = [custom_redblue; custom_red(idx_color(end:-1:1), :)];

list_ticklabel = cell(nLevel*2,1);
idx = 0;
for i = nLevel:-1:1
    idx = idx + 1;
    list_ticklabel{idx} = sprintf('p<%g', list_pval(i));
end
for i = 1:nLevel
    idx = idx + 1;
    list_ticklabel{idx} = sprintf('p<%g', list_pval(i));
end
list_tick = linspace(0,1,nLevel*2+2);
list_tick(end) = [];
list_tick = list_tick + (list_tick(2)-list_tick(1))/2;
list_tick(nLevel+1) = [];



%%%%% color bar %%%%%
figure;
fg = fig_setting_default();
fg.position = get(gcf, 'Position');
fg.position(3) = fg.position(3)*0.5;
set(gcf, 'Position', fg.position);
set(gcf, 'PaperPosition', fg.position);
axis off
custom_color = custom_redblue;
colormap(custom_color);
c = colorbar('fontsize',24, 'Ticks', list_tick, 'TickLabels', list_ticklabel);
c.Label.String = sprintf('Decrease              Increase');
c.LineWidth = 2;
c.Position(1) = 0.5;
c.Label.Position = [-2,0.5,0];
% save figure
outputfile = fullfile(dirFig, 'colorbar_pval_posneg');
print(outputfile,'-depsc','-opengl');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% task-event effects %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear lesion_result
list_factor = {'CPP', 'RU', 'Reward', 'Residual'};
nFactor = numel(list_factor);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(nFactor, nInteraction);
    lesion_result.(coef_name).sem = NaN(nFactor, nInteraction);
    lesion_result.(coef_name).pval = NaN(nFactor, nInteraction);
    lesion_result.(coef_name).effect = NaN(nFactor, nInteraction);
    
    % orignal data
    current_data.all = allData.all.reg.coef_behav.cpp_ru_vals.(coef_name);
    
    for idx_interaction = 1:nInteraction
        
        % lesion data
        current_data.lesion = allData.lesion{idx_interaction}.reg.coef_behav.cpp_ru_vals.(coef_name);
        
        x1 = current_data.lesion(:,:,idx_subgraph);
        x2 = current_data.all(:,:,idx_subgraph);
        
        data_diff = x1-x2; % between_lesion - between
        
        nSub = size(data_diff,1);
        
        diff_mean = mean(data_diff);
        diff_std = std(data_diff);
        diff_sem = diff_std./sqrt(nSub);
        
        [h,pval,ci,stats] = ttest(data_diff);
        
        lesion_result.(coef_name).raw(:, :, idx_interaction) = data_diff;
        lesion_result.(coef_name).mean(:,idx_interaction) = diff_mean;
        lesion_result.(coef_name).sem(:,idx_interaction) = diff_sem;
        lesion_result.(coef_name).pval(:,idx_interaction) = pval;
        lesion_result.(coef_name).effect(:,idx_interaction) = sign(stats.tstat);
        
    end
    
end


% graph
pval_threshold = 0.05/91;

coef_name = 'rel';
for idx_factor = 1:nFactor
    
    % factor_name
    factor_name = list_factor{idx_factor};
    
    mat_mean = zeros(nSystem, nSystem);
    mat_pval = zeros(nSystem, nSystem);
    idx_interaction = 0;
    for i = 1:nSystem
        for j = i:nSystem
            
            idx_interaction = idx_interaction + 1;
            mat_mean(i,j) = lesion_result.(coef_name).mean(idx_factor, idx_interaction);
            mat_pval(i,j) = lesion_result.(coef_name).pval(idx_factor, idx_interaction);
            
        end
    end
    mat_mean_original = mat_mean;
    mat_pval_original = mat_pval;
    
    
    %%%%% difference %%%%%
    mat_mean = mat_mean_original;
    mat_pval = mat_pval_original;
    
    idx_diag = logical(eye(nSystem));
    
    val_diag = mat_mean(idx_diag);
    mat_mean = mat_mean + mat_mean';
    mat_mean(idx_diag) = val_diag;
    
    val_diag = mat_pval(idx_diag);
    mat_pval = mat_pval + mat_pval';
    mat_pval(idx_diag) = val_diag;
    
    idx_sig = (mat_pval<pval_threshold);
    mat_mean(~idx_sig) = 0;
    mat_pval(~idx_sig) = 1;
    
    idx_pos = mat_mean>0;
    idx_neg = mat_mean<0;
    
    mat_mean_pos = mat_mean;
    mat_mean_neg = mat_mean;
    mat_mean_pos(idx_neg) = 0;
    mat_mean_neg(idx_pos) = 0;
    
    mat_pval_pos = mat_pval;
    mat_pval_neg = mat_pval;
    mat_pval_pos(idx_neg) = 1;
    mat_pval_neg(idx_pos) = 1;
    
    
    for i = 1:2
        
        switch i
            case 1
                switch factor_name
                    case {'CPP'}
                        clim = [0, 0.005];
                        ticks = [0:0.001:0.005];
                    case {'RU'}
                        clim = [0, 0.025];
                        ticks = [0:0.005:0.025];
                    case {'Reward'}
                        clim = [0, 0.005];
                        ticks = [0:0.001:0.005];
                    case {'Residual'}
                        clim = [0, 0.005];
                        ticks = [0:0.001:0.005];
                end
                mat_data = mat_mean_pos;
                custom_color = custom_parula;
                custom_color(1,:) = [1,1,1];
            case 2
                switch factor_name
                    case {'CPP'}
                        clim = [-0.005, 0];
                        ticks = [-0.005, -0.004, -0.003, -0.002, -0.001, 0];
                    case {'RU'}
                        clim = [-0.025, 0];
                        ticks = [-0.025, -0.02, -0.015, -0.01, -0.005, 0];
                    case {'Reward'}
                        clim = [-0.005, 0];
                        ticks = [-0.005, -0.004, -0.003, -0.002, -0.001, 0];
                    case {'Residual'}
                        clim = [-0.005, 0];
                        ticks = [-0.005, -0.004, -0.003, -0.002, -0.001, 0];
                end
                mat_data = mat_mean_neg;
                custom_color = custom_parula;
                custom_color(1,:) = [1,1,1];
                custom_color(1:end,:) = custom_color(end:-1:1,:);
        end
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*1.2;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        imagesc(mat_data, clim);
        
        colormap(custom_color);
        c = colorbar('fontsize',24, 'Ticks', ticks, 'TickLabels', ticks);
        c.Label.String = 'Change of coefficient';
        c.LineWidth = 2;
        c.Position(1) = 0.82;
        c.Label.Position(1) = -1.6;
        
        switch i
            case 1
                title(sprintf('%s', factor_name), 'fontsize',24);
                effect_type = 'pos';
            case 2
                c.Direction = 'reverse';
                title(sprintf('%s', factor_name), 'fontsize',24);
                effect_type = 'neg';
        end
        
        set(gca,...
            'XTick',[1:2:13],...
            'YTick',[1:2:13],...
            'fontsize',24,...
            'linewidth',2,...
            'box','on');
        
        xlabel('Systems', 'fontsize', 24);
        ylabel('Systems', 'fontsize', 24);
        
        axis image
        
        % outputfile
        outputfile = fullfile(dirFig, sprintf('coef_effect_lesion_interaction_%s_subgraph%.2d_%s_%s',coef_name, idx_subgraph, factor_name, effect_type));
        print(outputfile,'-depsc','-opengl');
        
    end
    
    
    
    %%%%% pval %%%%%
    mat_mean = mat_mean_original;
    mat_pval = mat_pval_original;
    
    idx_diag = logical(eye(nSystem));
    
    val_diag = mat_mean(idx_diag);
    mat_mean = mat_mean + mat_mean';
    mat_mean(idx_diag) = val_diag;
    
    val_diag = mat_pval(idx_diag);
    mat_pval = mat_pval + mat_pval';
    mat_pval(idx_diag) = val_diag;
    
    
    mat_pval_idx = NaN(size(nSystem, nSystem));
    for i = 1:nSystem
        for j = 1:nSystem
            
            val = mat_mean(i,j);
            pval = mat_pval(i,j);
            if val>=0
                if pval<list_pval(1) & pval>=list_pval(2)
                    pval_idx = 6;
                elseif pval<list_pval(2) & pval>=list_pval(3)
                    pval_idx = 7;
                elseif pval<list_pval(3) & pval>=list_pval(4)
                    pval_idx = 8;
                elseif pval<list_pval(4)
                    pval_idx = 9;
                else
                    pval_idx = 5;
                end
            else
                if pval<list_pval(1) & pval>=list_pval(2)
                    pval_idx = 4;
                elseif pval<list_pval(2) & pval>=list_pval(3)
                    pval_idx = 3;
                elseif pval<list_pval(3) & pval>=list_pval(4)
                    pval_idx = 2;
                elseif pval<list_pval(4)
                    pval_idx = 1;
                else
                    pval_idx = 5;
                end
            end
            
            mat_pval_idx(i,j) = pval_idx;
            
        end
    end
    
    unique_idx = unique(mat_pval_idx(:));
    custom_color = custom_redblue(unique_idx,:);
    mat_data = mat_pval_idx;
    for i = 1:numel(unique_idx)
        idx_replace = (mat_pval_idx==unique_idx(i));
        mat_pval_idx(idx_replace) = i;
    end
    
    figure;
    fg = fig_setting_default;
    
    imagesc(mat_pval_idx);
    colormap(custom_color);
    
    title(sprintf('%s', factor_name), 'fontsize',24);
    
    set(gca,...
        'XTick',[1:2:13],...
        'YTick',[1:2:13],...
        'fontsize',24,...
        'linewidth',2,...
        'box','on');
    
    xlabel('Systems', 'fontsize', 24);
    ylabel('Systems', 'fontsize', 24);
    
    axis image
    
    % highlight significant subgraph system by system edges
    hold on
    for i = 1:nSystem
        for j = 1:nSystem
            
            if subnet_thresholded(i,j)==1
                
                plot(i,j,...
                    'color', 'k',...
                    'linestyle', 'none',...
                    'linewidth', 2,...
                    'marker', 'o',...
                    'markersize', 20);
                
            end
            
        end
    end
    hold off
    
    % outputfile
    outputfile = fullfile(dirFig, sprintf('coef_effect_lesion_interaction_%s_subgraph%.2d_%s_pval_posneg',coef_name, idx_subgraph, factor_name));
    print(outputfile,'-depsc','-opengl');
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% individual difference %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normative learning and dynamic modulation
tic
clear lesion_result;
normative_LR = allData.all.behav.behav_beta(:,3) + allData.all.behav.behav_beta(:,4);
list_factor = {'Dynamic modulation'};
nFactor = numel(list_factor);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(1, nInteraction);
    lesion_result.(coef_name).sem = NaN(1, nInteraction);
    lesion_result.(coef_name).pval = NaN(1, nInteraction);
    
    % orignal data
    current_data.all = allData.all.reg.coef_behav.cpp_ru_vals.(coef_name)(:,1,:) +...
        allData.all.reg.coef_behav.cpp_ru_vals.(coef_name)(:,2,:);
    
    for idx_interaction = 1:nInteraction
        
        % lesion data
        current_data.lesion = allData.lesion{idx_interaction}.reg.coef_behav.cpp_ru_vals.(coef_name)(:,1,:) +...
            allData.lesion{idx_interaction}.reg.coef_behav.cpp_ru_vals.(coef_name)(:,2,:);
        
        x1 = current_data.lesion(:,:,idx_subgraph);
        x2 = current_data.all(:,:,idx_subgraph);
        
        
        [rho_x1] = corr(normative_LR, x1);
        [rho_x2] = corr(normative_LR, x2);
        
        diff_mean = rho_x1 - rho_x2;
        
        nSub = size(x1,1);
        temp_diff = NaN(nPermutation,1);
        for p = 1:nPermutation
            
            idx_random = randperm(nSub);
            normative_LR_perm = normative_LR(idx_random);
            [rho_x1_perm] = corr(normative_LR_perm, x1);
            [rho_x2_perm] = corr(normative_LR_perm, x2);
            diff_mean_perm = rho_x1_perm - rho_x2_perm;
            temp_diff(p,1) = diff_mean_perm;
            
        end
        diff_sem = std(temp_diff);
        
        if diff_mean>=0
            pval = mean(temp_diff>diff_mean);
        else
            pval = mean(temp_diff<diff_mean);
        end
        
        
        lesion_result.(coef_name).mean(1,idx_interaction) = diff_mean;
        lesion_result.(coef_name).sem(1,idx_interaction) = diff_sem;
        lesion_result.(coef_name).pval(1,idx_interaction) = pval;
        
    end
    
end
toc




% graph
pval_threshold = 0.05/91;
coef_name = 'rel';
for idx_factor = 1:nFactor
    
    % factor_name
    factor_name = list_factor{idx_factor};
    
    mat_mean = zeros(nSystem, nSystem);
    mat_pval = zeros(nSystem, nSystem);
    idx_interaction = 0;
    for i = 1:nSystem
        for j = i:nSystem
            
            idx_interaction = idx_interaction + 1;
            mat_mean(i,j) = lesion_result.(coef_name).mean(idx_factor, idx_interaction);
            mat_pval(i,j) = lesion_result.(coef_name).pval(idx_factor, idx_interaction);
            
        end
    end
    mat_mean_original = mat_mean;
    mat_pval_original = mat_pval;
    
    
    %%%%% difference %%%%%
    mat_mean = mat_mean_original;
    mat_pval = mat_pval_original;
    
    idx_diag = logical(eye(nSystem));
    
    val_diag = mat_mean(idx_diag);
    mat_mean = mat_mean + mat_mean';
    mat_mean(idx_diag) = val_diag;
    
    val_diag = mat_pval(idx_diag);
    mat_pval = mat_pval + mat_pval';
    mat_pval(idx_diag) = val_diag;
    
    idx_sig = (mat_pval<pval_threshold);
    mat_mean(~idx_sig) = 0;
    mat_pval(~idx_sig) = 1;
    
    idx_pos = mat_mean>0;
    idx_neg = mat_mean<0;
    
    mat_mean_pos = mat_mean;
    mat_mean_neg = mat_mean;
    mat_mean_pos(idx_neg) = 0;
    mat_mean_neg(idx_pos) = 0;
    
    mat_pval_pos = mat_pval;
    mat_pval_neg = mat_pval;
    mat_pval_pos(idx_neg) = 1;
    mat_pval_neg(idx_pos) = 1;
    
    for i = 1:2
        
        switch i
            case 1
                clim = [0, 0.06];
                ticks = [0:0.02:0.06];
                mat_data = mat_mean_pos;
                custom_color = custom_parula;
                custom_color(1,:) = [1,1,1];
            case 2
                clim = [-0.04, 0];
                ticks = [-0.04:0.01:0];
                mat_data = mat_mean_neg;
                custom_color = custom_parula;
                custom_color(1,:) = [1,1,1];
                custom_color(1:end,:) = custom_color(end:-1:1,:);
        end
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*1.2;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        imagesc(mat_data, clim);
        
        colormap(custom_color);
        c = colorbar('fontsize',24, 'Ticks', ticks, 'TickLabels', ticks);
        c.Label.String = 'Change of coefficient';
        c.LineWidth = 2;
        c.Position(1) = 0.82;
        c.Label.Position(1) = -1.6;
        switch i
            case 1
                title(sprintf('%s', factor_name), 'fontsize',24);
                effect_type = 'pos';
            case 2
                c.Direction = 'reverse';
                title(sprintf('%s', factor_name), 'fontsize',24);
                effect_type = 'neg';
        end
        
        set(gca,...
            'XTick',[1:2:13],...
            'YTick',[1:2:13],...
            'fontsize',24,...
            'linewidth',2,...
            'box','on');
        
        xlabel('Systems', 'fontsize', 24);
        ylabel('Systems', 'fontsize', 24);
        
        axis image
        
        % outputfile
        outputfile = fullfile(dirFig, sprintf('coef_corr_dynamic_lesion_interaction_%s_subgraph%.2d_%s',coef_name, idx_subgraph, effect_type));
        print(outputfile,'-depsc','-opengl');
        
        
    end
    
    %%%%% pval %%%%%
    mat_mean = mat_mean_original;
    mat_pval = mat_pval_original;
    
    idx_diag = logical(eye(nSystem));
    
    val_diag = mat_mean(idx_diag);
    mat_mean = mat_mean + mat_mean';
    mat_mean(idx_diag) = val_diag;
    
    val_diag = mat_pval(idx_diag);
    mat_pval = mat_pval + mat_pval';
    mat_pval(idx_diag) = val_diag;
    
    
    mat_pval_idx = NaN(size(nSystem, nSystem));
    for i = 1:nSystem
        for j = 1:nSystem
            
            val = mat_mean(i,j);
            pval = mat_pval(i,j);
            if val>=0
                if pval<list_pval(1) & pval>=list_pval(2)
                    pval_idx = 6;
                elseif pval<list_pval(2) & pval>=list_pval(3)
                    pval_idx = 7;
                elseif pval<list_pval(3) & pval>=list_pval(4)
                    pval_idx = 8;
                elseif pval<list_pval(4)
                    pval_idx = 9;
                else
                    pval_idx = 5;
                end
            else
                if pval<list_pval(1) & pval>=list_pval(2)
                    pval_idx = 4;
                elseif pval<list_pval(2) & pval>=list_pval(3)
                    pval_idx = 3;
                elseif pval<list_pval(3) & pval>=list_pval(4)
                    pval_idx = 2;
                elseif pval<list_pval(4)
                    pval_idx = 1;
                else
                    pval_idx = 5;
                end
            end
            
            mat_pval_idx(i,j) = pval_idx;
            
        end
    end
    
    unique_idx = unique(mat_pval_idx(:));
    custom_color = custom_redblue(unique_idx,:);
    mat_data = mat_pval_idx;
    for i = 1:numel(unique_idx)
        idx_replace = (mat_pval_idx==unique_idx(i));
        mat_pval_idx(idx_replace) = i;
    end
    
    figure;
    fg = fig_setting_default;
    
    imagesc(mat_pval_idx);
    colormap(custom_color);
    
    title(sprintf('%s', factor_name), 'fontsize',24);
    
    set(gca,...
        'XTick',[1:2:13],...
        'YTick',[1:2:13],...
        'fontsize',24,...
        'linewidth',2,...
        'box','on');
    
    xlabel('Systems', 'fontsize', 24);
    ylabel('Systems', 'fontsize', 24);
    
    axis image
    
    highlight significant subgraph system by system edges
    hold on
    for i = 1:nSystem
        for j = 1:nSystem
            
            if subnet_thresholded(i,j)==1
                
                plot(i,j,...
                    'color', 'k',...
                    'linestyle', 'none',...
                    'linewidth', 2,...
                    'marker', 'o',...
                    'markersize', 20);
                
            end
            
        end
    end
    hold off
    
    % outputfile
    outputfile = fullfile(dirFig, sprintf('coef_corr_dynamic_lesion_interaction_%s_subgraph%.2d_pval_posneg',coef_name, idx_subgraph));
    print(outputfile,'-depsc','-opengl');
    
    
    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% individual difference %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normative learning and average expression
tic
clear lesion_result;
normative_LR = allData.all.behav.behav_beta(:,3) + allData.all.behav.behav_beta(:,4);
list_factor = {'Average expression'};
nFactor = numel(list_factor);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(1, nInteraction);
    lesion_result.(coef_name).sem = NaN(1, nInteraction);
    lesion_result.(coef_name).pval = NaN(1, nInteraction);
    
    % orignal data
    current_data.all = allData.all.subgraph.subject.(coef_name).mean;
    
    for idx_interaction = 1:nInteraction
        
        % lesion data
        current_data.lesion = allData.lesion{idx_interaction}.subgraph.subject.(coef_name).mean;
        
        x1 = current_data.lesion(:,idx_subgraph);
        x2 = current_data.all(:,idx_subgraph);
        
        
        [rho_x1] = corr(normative_LR, x1);
        [rho_x2] = corr(normative_LR, x2);
        
        diff_mean = rho_x1 - rho_x2;
        
        nSub = size(x1,1);
        temp_diff = NaN(nPermutation,1);
        for p = 1:nPermutation
            
            idx_random = randperm(nSub);
            normative_LR_perm = normative_LR(idx_random);
            [rho_x1_perm] = corr(normative_LR_perm, x1);
            [rho_x2_perm] = corr(normative_LR_perm, x2);
            diff_mean_perm = rho_x1_perm - rho_x2_perm;
            temp_diff(p,1) = diff_mean_perm;
            
        end
        diff_sem = std(temp_diff);
        
        if diff_mean>=0
            pval = mean(temp_diff>diff_mean);
        else
            pval = mean(temp_diff<diff_mean);
        end
        
        
        lesion_result.(coef_name).mean(1,idx_interaction) = diff_mean;
        lesion_result.(coef_name).sem(1,idx_interaction) = diff_sem;
        lesion_result.(coef_name).pval(1,idx_interaction) = pval;
        
    end
    
end
toc



% graph
pval_threshold = 0.05/91;
coef_name = 'rel';
for idx_factor = 1:nFactor
    
    % factor_name
    factor_name = list_factor{idx_factor};
    
    mat_mean = zeros(nSystem, nSystem);
    mat_pval = zeros(nSystem, nSystem);
    idx_interaction = 0;
    for i = 1:nSystem
        for j = i:nSystem
            
            idx_interaction = idx_interaction + 1;
            mat_mean(i,j) = lesion_result.(coef_name).mean(idx_factor, idx_interaction);
            mat_pval(i,j) = lesion_result.(coef_name).pval(idx_factor, idx_interaction);
            
        end
    end
    mat_mean_original = mat_mean;
    mat_pval_original = mat_pval;
    
    
    %%%%% difference %%%%%
    mat_mean = mat_mean_original;
    mat_pval = mat_pval_original;
    
    idx_diag = logical(eye(nSystem));
    
    val_diag = mat_mean(idx_diag);
    mat_mean = mat_mean + mat_mean';
    mat_mean(idx_diag) = val_diag;
    
    val_diag = mat_pval(idx_diag);
    mat_pval = mat_pval + mat_pval';
    mat_pval(idx_diag) = val_diag;
    
    idx_sig = (mat_pval<pval_threshold);
    mat_mean(~idx_sig) = 0;
    mat_pval(~idx_sig) = 1;
    
    idx_pos = mat_mean>0;
    idx_neg = mat_mean<0;
    
    mat_mean_pos = mat_mean;
    mat_mean_neg = mat_mean;
    mat_mean_pos(idx_neg) = 0;
    mat_mean_neg(idx_pos) = 0;
    
    mat_pval_pos = mat_pval;
    mat_pval_neg = mat_pval;
    mat_pval_pos(idx_neg) = 1;
    mat_pval_neg(idx_pos) = 1;
    
    
    for i = 1:2
        
        switch i
            case 1
                clim = [0, 0.008];
                ticks = [0:0.002:0.008];
                mat_data = mat_mean_pos;
                custom_color = custom_parula;
                custom_color(1,:) = [1,1,1];
            case 2
                clim = [-0.008, 0];
                ticks = [-0.008:0.002:0];
                mat_data = mat_mean_neg;
                custom_color = custom_parula;
                custom_color(1,:) = [1,1,1];
                custom_color(1:end,:) = custom_color(end:-1:1,:);
        end
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*1.2;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        imagesc(mat_data, clim);
        
        colormap(custom_color);
        c = colorbar('fontsize',24, 'Ticks', ticks, 'TickLabels', ticks);
        c.Label.String = 'Change of coefficient';
        c.LineWidth = 2;
        c.Position(1) = 0.82;
        c.Label.Position(1) = -1.6;
        switch i
            case 1
                title(sprintf('%s', factor_name), 'fontsize',24);
                effect_type = 'pos';
            case 2
                c.Direction = 'reverse';
                title(sprintf('%s', factor_name), 'fontsize',24);
                effect_type = 'neg';
        end
        
        set(gca,...
            'XTick',[1:2:13],...
            'YTick',[1:2:13],...
            'fontsize',24,...
            'linewidth',2,...
            'box','on');
        
        xlabel('Systems', 'fontsize', 24);
        ylabel('Systems', 'fontsize', 24);
        
        axis image
        
        
        % outputfile
        outputfile = fullfile(dirFig, sprintf('coef_corr_average_lesion_interaction_%s_subgraph%.2d_%s',coef_name, idx_subgraph, effect_type));
        print(outputfile,'-depsc','-opengl');
        
    end
    
    
    
    %%%%% pval %%%%%
    mat_mean = mat_mean_original;
    mat_pval = mat_pval_original;
    
    idx_diag = logical(eye(nSystem));
    
    val_diag = mat_mean(idx_diag);
    mat_mean = mat_mean + mat_mean';
    mat_mean(idx_diag) = val_diag;
    
    val_diag = mat_pval(idx_diag);
    mat_pval = mat_pval + mat_pval';
    mat_pval(idx_diag) = val_diag;
    
    
    mat_pval_idx = NaN(size(nSystem, nSystem));
    for i = 1:nSystem
        for j = 1:nSystem
            
            val = mat_mean(i,j);
            pval = mat_pval(i,j);
            if val>=0
                if pval<list_pval(1) & pval>=list_pval(2)
                    pval_idx = 6;
                elseif pval<list_pval(2) & pval>=list_pval(3)
                    pval_idx = 7;
                elseif pval<list_pval(3) & pval>=list_pval(4)
                    pval_idx = 8;
                elseif pval<list_pval(4)
                    pval_idx = 9;
                else
                    pval_idx = 5;
                end
            else
                if pval<list_pval(1) & pval>=list_pval(2)
                    pval_idx = 4;
                elseif pval<list_pval(2) & pval>=list_pval(3)
                    pval_idx = 3;
                elseif pval<list_pval(3) & pval>=list_pval(4)
                    pval_idx = 2;
                elseif pval<list_pval(4)
                    pval_idx = 1;
                else
                    pval_idx = 5;
                end
            end
            
            mat_pval_idx(i,j) = pval_idx;
            
        end
    end
    
    unique_idx = unique(mat_pval_idx(:));
    custom_color = custom_redblue(unique_idx,:);
    mat_data = mat_pval_idx;
    for i = 1:numel(unique_idx)
        idx_replace = (mat_pval_idx==unique_idx(i));
        mat_pval_idx(idx_replace) = i;
    end
    
    figure;
    fg = fig_setting_default;
    
    imagesc(mat_pval_idx);
    colormap(custom_color);
    
    title(sprintf('%s', factor_name), 'fontsize',24);
    
    set(gca,...
        'XTick',[1:2:13],...
        'YTick',[1:2:13],...
        'fontsize',24,...
        'linewidth',2,...
        'box','on');
    
    xlabel('Systems', 'fontsize', 24);
    ylabel('Systems', 'fontsize', 24);
    
    axis image
    
    % highlight significant subgraph system by system edges
    hold on
    for i = 1:nSystem
        for j = 1:nSystem
            
            if subnet_thresholded(i,j)==1
                
                plot(i,j,...
                    'color', 'k',...
                    'linestyle', 'none',...
                    'linewidth', 2,...
                    'marker', 'o',...
                    'markersize', 20);
                
            end
            
        end
    end
    hold off
    
    % outputfile
    outputfile = fullfile(dirFig, sprintf('coef_corr_average_lesion_interaction_%s_subgraph%.2d_pval_posneg',coef_name, idx_subgraph));
    print(outputfile,'-depsc','-opengl');
    
    
    
end



