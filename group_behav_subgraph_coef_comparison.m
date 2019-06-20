function group_behav_subgraph_coef_comparison()

window_width = 10;
window_overlap = 8;

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirFig = fullfile(dirType, sprintf('figures_%d_%d', window_width, window_overlap));
mkdir(dirFig);

% setting
nPermutation = 10000;
list_coef = {'rel'};
nCoef = numel(list_coef);


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
    'Dorsal attention';
    };
nSystem = numel(list_system);


% load data
allData.all = group_behav_subgraph_coef('', 0);
allData.within = group_behav_subgraph_coef('_within', 0);
allData.between = group_behav_subgraph_coef('_between', 0);
for idx_system = 1:nSystem
    allData.lesion{idx_system} = group_behav_subgraph_coef(sprintf('_lesion_%d',idx_system),0);
end
for idx_system = 1:nSystem
    allData.between_lesion{idx_system} = group_behav_subgraph_coef(sprintf('_between_lesion_%d',idx_system),0);
end


nSub = size(allData.all.reg.coef_behav.cpp_ru_vals.rel,1);

idx_subgraph = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% system: lesion %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% task-event effects
clear lesion_result
list_factor = {'CPP', 'RU', 'Reward', 'Residual'};
nFactor = numel(list_factor);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(nFactor, nSystem);
    lesion_result.(coef_name).sem = NaN(nFactor, nSystem);
    lesion_result.(coef_name).pval = NaN(nFactor, nSystem);
    lesion_result.(coef_name).effect = NaN(nFactor, nSystem);
    
    % orignal data
    current_data.all = allData.all.reg.coef_behav.cpp_ru_vals.(coef_name);
    
    for idx_system = 1:nSystem
        
        % lesion data
        current_data.lesion = allData.lesion{idx_system}.reg.coef_behav.cpp_ru_vals.(coef_name);
        
        x1 = current_data.lesion(:,:,idx_subgraph);
        x2 = current_data.all(:,:,idx_subgraph);
        
        data_diff = x1-x2; % lesion - all
        
        nSub = size(data_diff,1);
        
        diff_mean = mean(data_diff);
        diff_std = std(data_diff);
        diff_sem = diff_std./sqrt(nSub);
        
        [h,pval,ci,stats] = ttest(data_diff);
        
        lesion_result.(coef_name).mean(:,idx_system) = diff_mean;
        lesion_result.(coef_name).sem(:,idx_system) = diff_sem;
        lesion_result.(coef_name).pval(:,idx_system) = pval;
        lesion_result.(coef_name).effect(:,idx_system) = sign(stats.tstat);
        
    end
    
end

plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
plotParameter.color.bar.edgecolor = {'none'};
plotParameter.color.error = [0, 0, 0];
plotParameter.linewidth = 2;


threshold_pval = .05;
for c = 1:nCoef
    
    coef_name = list_coef{c};
    
    current_mean = lesion_result.(coef_name).mean;
    current_sem = lesion_result.(coef_name).sem;
    current_pval = lesion_result.(coef_name).pval;
    current_effect = lesion_result.(coef_name).effect;
    current_h = double(current_pval<threshold_pval);
    
    %%%%% horizontal bar %%%%%
    for f = 1:nFactor
        
        factor_name = list_factor{f};
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*0.9;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        
        meanData = current_mean(f,:);
        semData = current_sem(f,:);
        pvalData = current_pval(f,:);
        
        hold on
        
        plot([0;0], [0.5;numel(meanData)+0.5],...
            'linestyle', '--',...
            'linewidth', 2,...
            'color', [0.5,0.5,0.5]);
        e = errorbar(meanData', [1:nSystem]', semData', 'horizontal');
        set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
        plot(meanData', [1:nSystem]',...
            'marker', 'o',...
            'linestyle', 'none',...
            'linewidth', 2,...
            'markersize', 12,...
            'markerfacecolor', [0.5,0.5,0.5],...
            'markeredgecolor', 'k');
        
        hold off
        
        set(gca,'YTick',[1:nSystem],'YTickLabel',list_system);
        set(gca, 'fontsize',18, 'linewidth',2);
        
        xlabel('Change of coefficient','fontsize',24);
        ylabel('Lesion system','fontsize',24);
        
        ylim([0.5, numel(meanData)+0.5]);
        yrange = ylim;
        
        switch factor_name
            case {'CPP'}
                xrange = [-0.08, 0.08];
            case {'RU'}
                xrange = [-0.15, 0.15];
            case {'Reward'}
                xrange = [-0.04, 0.04];
            case {'Residual'}
                xrange = [-0.08, 0.08];
        end
        xlim(xrange);
        
        
        for i = 1:numel(pvalData)
            
            idx_p = 0;
            if pvalData(i)<.05 & pvalData(i)>=0.01
                idx_p = 1;
                list_n = 0;
            elseif pvalData(i)<0.01 & pvalData(i)>=0.001
                idx_p = 1;
                list_n = [-0.5,0.5];
            elseif pvalData(i)<0.001
                idx_p = 1;
                list_n = [-1,0,1];
            end
            
            if idx_p==1
                y = i;
                if meanData(i)>=0
                    x = meanData(i)+semData(i) + 0.05*(max(xrange)-min(xrange));
                else
                    x = meanData(i)-semData(i) - 0.05*(max(xrange)-min(xrange));
                end
                
                list_n = list_n*0.27;
                for n = 1:numel(list_n)
                    newy = y + list_n(n) - 0.2;
                    text(x,newy,'*','HorizontalAlignment','left','VerticalAlignment', 'middle', 'fontsize',24);
                end
                
            end
            
        end
        
        title(factor_name, 'fontsize', 24);
        
        outputfile = fullfile(dirFig, sprintf('coef_effect_lesion_system_%s_subgraph%.2d_%s',coef_name, idx_subgraph, factor_name));
        print(outputfile,'-depsc','-opengl');
        
    end
    
end


% individual difference of normative learning and dynamic modulation
clear lesion_result;
normative_LR = allData.all.behav.behav_beta(:,3) + allData.all.behav.behav_beta(:,4);
list_factor = {'Dynamic modulation'};
nFactor = numel(list_factor);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(1, nSystem);
    lesion_result.(coef_name).sem = NaN(1, nSystem);
    lesion_result.(coef_name).pval = NaN(1, nSystem);
    
    % orignal data
    current_data.all = allData.all.reg.coef_behav.cpp_ru_vals.(coef_name)(:,1,:) +...
        allData.all.reg.coef_behav.cpp_ru_vals.(coef_name)(:,2,:);
    
    for idx_system = 1:nSystem
        
        % lesion data
        current_data.lesion = allData.lesion{idx_system}.reg.coef_behav.cpp_ru_vals.(coef_name)(:,1,:) +...
            allData.lesion{idx_system}.reg.coef_behav.cpp_ru_vals.(coef_name)(:,2,:);
        
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
        
        
        lesion_result.(coef_name).mean(1,idx_system) = diff_mean;
        lesion_result.(coef_name).sem(1,idx_system) = diff_sem;
        lesion_result.(coef_name).pval(1,idx_system) = pval;
        
    end
    
end

plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
plotParameter.color.bar.edgecolor = {'none'};
plotParameter.color.error = [0, 0, 0];
plotParameter.linewidth = 2;

threshold_pval = .05;
for c = 1:nCoef
    
    coef_name = list_coef{c};
    
    current_mean = lesion_result.(coef_name).mean(1,:);
    current_sem = lesion_result.(coef_name).sem(1,:);
    current_pval = lesion_result.(coef_name).pval(1,:);
    
    %%%%% horizontal bar %%%%%
    for f = 1:nFactor
        
        factor_name = list_factor{f};
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*0.9;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        
        meanData = current_mean(f,:);
        semData = current_sem(f,:);
        pvalData = current_pval(f,:);
        
        hold on
        
        plot([0;0], [0.5;numel(meanData)+0.5],...
            'linestyle', '--',...
            'linewidth', 2,...
            'color', [0.5,0.5,0.5]);
        e = errorbar(meanData', [1:nSystem]', semData', 'horizontal');
        set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
        plot(meanData', [1:nSystem]',...
            'marker', 'o',...
            'linestyle', 'none',...
            'linewidth', 2,...
            'markersize', 12,...
            'markerfacecolor', [0.5,0.5,0.5],...
            'markeredgecolor', 'k');
        
        hold off
        
        set(gca,'YTick',[1:nSystem],'YTickLabel',list_system);
        set(gca, 'fontsize',18, 'linewidth',2);
        
        xlabel('Change of coefficient','fontsize',24);
        ylabel('Lesion system','fontsize',24);
        
        ylim([0.5, numel(meanData)+0.5]);
        yrange = ylim;
        
        
        xrange = [-0.2, 0.2];
        xlim(xrange);
        
        
        for i = 1:numel(pvalData)
            
            idx_p = 0;
            if pvalData(i)<.05 & pvalData(i)>=0.01
                idx_p = 1;
                list_n = 0;
            elseif pvalData(i)<0.01 & pvalData(i)>=0.001
                idx_p = 1;
                list_n = [-0.5,0.5];
            elseif pvalData(i)<0.001
                idx_p = 1;
                list_n = [-1,0,1];
            end
            
            if idx_p==1
                y = i;
                if meanData(i)>=0
                    x = meanData(i)+semData(i) + 0.05*(max(xrange)-min(xrange));
                else
                    x = meanData(i)-semData(i) - 0.05*(max(xrange)-min(xrange));
                end
                
                list_n = list_n*0.27;
                for n = 1:numel(list_n)
                    newy = y + list_n(n) - 0.2;
                    text(x,newy,'*','HorizontalAlignment','left','VerticalAlignment', 'middle', 'fontsize',24);
                end
                
            end
            
        end
        
        title(factor_name, 'fontsize', 24);
        
        outputfile = fullfile(dirFig, sprintf('coef_corr_dynamic_lesion_system_%s_subgraph%.2d',coef_name, idx_subgraph));
        print(outputfile,'-depsc','-opengl');
        
    end
    
end


% individual difference of normative learning
clear lesion_result;
list_factor = {'Average expression'};
nFactor = numel(list_factor);
normative_LR = allData.all.behav.behav_beta(:,3) + allData.all.behav.behav_beta(:,4);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(1, nSystem);
    lesion_result.(coef_name).sem = NaN(1, nSystem);
    lesion_result.(coef_name).pval = NaN(1, nSystem);
    
    % orignal data
    current_data.all = allData.all.subgraph.subject.(coef_name).mean;
    
    for idx_system = 1:nSystem
        
        % lesion data
        current_data.lesion = allData.lesion{idx_system}.subgraph.subject.(coef_name).mean;
        
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
        
        
        lesion_result.(coef_name).mean(1,idx_system) = diff_mean;
        lesion_result.(coef_name).sem(1,idx_system) = diff_sem;
        lesion_result.(coef_name).pval(1,idx_system) = pval;
        
    end
    
end

plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
plotParameter.color.bar.edgecolor = {'none'};
plotParameter.color.error = [0, 0, 0];
plotParameter.linewidth = 2;

threshold_pval = .05;
for c = 1:nCoef
    
    coef_name = list_coef{c};
    
    current_mean = lesion_result.(coef_name).mean(1,:);
    current_sem = lesion_result.(coef_name).sem(1,:);
    current_pval = lesion_result.(coef_name).pval(1,:);
    
    %%%%% horizontal bar %%%%%
    for f = 1:nFactor
        
        factor_name = list_factor{f};
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*0.9;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        
        meanData = current_mean(f,:);
        semData = current_sem(f,:);
        pvalData = current_pval(f,:);
        
        hold on
        
        plot([0;0], [0.5;numel(meanData)+0.5],...
            'linestyle', '--',...
            'linewidth', 2,...
            'color', [0.5,0.5,0.5]);
        e = errorbar(meanData', [1:nSystem]', semData', 'horizontal');
        set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
        plot(meanData', [1:nSystem]',...
            'marker', 'o',...
            'linestyle', 'none',...
            'linewidth', 2,...
            'markersize', 12,...
            'markerfacecolor', [0.5,0.5,0.5],...
            'markeredgecolor', 'k');
        
        hold off
        
        set(gca,'YTick',[1:nSystem],'YTickLabel',list_system);
        set(gca, 'fontsize',18, 'linewidth',2);
        
        xlabel('Change of coefficient','fontsize',24);
        ylabel('Lesion system','fontsize',24);
        
        ylim([0.5, numel(meanData)+0.5]);
        yrange = ylim;
        
        xrange = [-0.06, 0.06];
        xlim(xrange);
        
        
        for i = 1:numel(pvalData)
            
            idx_p = 0;
            if pvalData(i)<.05 & pvalData(i)>=0.01
                idx_p = 1;
                list_n = 0;
            elseif pvalData(i)<0.01 & pvalData(i)>=0.001
                idx_p = 1;
                list_n = [-0.5,0.5];
            elseif pvalData(i)<0.001
                idx_p = 1;
                list_n = [-1,0,1];
            end
            
            if idx_p==1
                y = i;
                if meanData(i)>=0
                    x = meanData(i)+semData(i) + 0.05*(max(xrange)-min(xrange));
                else
                    x = meanData(i)-semData(i) - 0.05*(max(xrange)-min(xrange));
                end
                
                list_n = list_n*0.27;
                for n = 1:numel(list_n)
                    newy = y + list_n(n) - 0.2;
                    text(x,newy,'*','HorizontalAlignment','left','VerticalAlignment', 'middle', 'fontsize',24);
                end
                
            end
            
        end
        
        title(factor_name, 'fontsize', 24);
        
        outputfile = fullfile(dirFig, sprintf('coef_corr_average_lesion_system_%s_subgraph%.2d',coef_name, idx_subgraph));
        print(outputfile,'-depsc','-opengl');
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% comparison: within versus between %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear lesion_result
list_comparison = {
    'within', 'all';
    'between', 'all';
    'between', 'within';
    };
nComparison = size(list_comparison,1);
list_type = {
    'Within vs All';
    'Between vs All';
    'Between vs Within';
    };

% tesk-event effects
list_factor = {'CPP', 'RU', 'Reward', 'Residual'};
nFactor = numel(list_factor);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(nFactor, nComparison);
    lesion_result.(coef_name).sem = NaN(nFactor, nComparison);
    lesion_result.(coef_name).pval = NaN(nFactor, nComparison);
    lesion_result.(coef_name).tstat = NaN(nFactor, nComparison);
    lesion_result.(coef_name).effect = NaN(nFactor, nComparison);
    
    % orignal data
    current_data.all = allData.all.reg.coef_behav.cpp_ru_vals.(coef_name);
    current_data.within = allData.within.reg.coef_behav.cpp_ru_vals.(coef_name);
    current_data.between = allData.between.reg.coef_behav.cpp_ru_vals.(coef_name);
    
    for idx_comparison = 1:nComparison
        
        x1 = current_data.(list_comparison{idx_comparison,1})(:,:,idx_subgraph);
        x2 = current_data.(list_comparison{idx_comparison,2})(:,:,idx_subgraph);
        
        data_diff = x1-x2;
        
        nSub = size(data_diff,1);
        
        diff_mean = mean(data_diff);
        diff_std = std(data_diff);
        diff_sem = diff_std./sqrt(nSub);
        
        [h,pval,ci,stats] = ttest(data_diff);
        
        lesion_result.(coef_name).raw(:,:,idx_comparison) = data_diff;
        lesion_result.(coef_name).mean(:,idx_comparison) = diff_mean;
        lesion_result.(coef_name).sem(:,idx_comparison) = diff_sem;
        lesion_result.(coef_name).pval(:,idx_comparison) = pval;
        lesion_result.(coef_name).tstat(:,idx_comparison) = stats.tstat;
        lesion_result.(coef_name).effect(:,idx_comparison) = sign(stats.tstat);
        
    end
    
end

% graph
plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
plotParameter.color.bar.edgecolor = {'none'};
plotParameter.color.error = [0, 0, 0];
plotParameter.linewidth = 2;

plotParameter.xticklabel = {'Within-All', 'Between-All', 'Within-Between'};

threshold_pval = .05;
for c = 1:nCoef
    
    coef_name = list_coef{c};
    
    current_mean = lesion_result.(coef_name).mean;
    current_sem = lesion_result.(coef_name).sem;
    current_pval = lesion_result.(coef_name).pval;
    current_effect = lesion_result.(coef_name).effect;
    current_h = double(current_pval<threshold_pval);
    current_disp = current_h.*current_effect;
    
    
    %%%%% horizontal bar %%%%%
    for f = 1:nFactor
        
        factor_name = list_factor{f};
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*0.85;
        fg.position(4) = fg.position(4)*0.45;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        
        meanData = current_mean(f,:);
        semData = current_sem(f,:);
        pvalData = current_pval(f,:);
        
        hold on
        
        plot([0;0], [0.5;numel(meanData)+0.5],...
            'linestyle', '--',...
            'linewidth', 2,...
            'color', [0.5,0.5,0.5]);
        e = errorbar(meanData', [1:nComparison]', semData', 'horizontal');
        set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
        plot(meanData', [1:nComparison]',...
            'marker', 'o',...
            'linestyle', 'none',...
            'linewidth', 2,...
            'markersize', 12,...
            'markerfacecolor', [0.5,0.5,0.5],...
            'markeredgecolor', 'k');
        
        hold off
        
        set(gca,'YTick',[1:nComparison],'YTickLabel',list_type);
        set(gca, 'fontsize',18, 'linewidth',2);
        
        xlabel('Change of coefficient','fontsize',24);
        ylabel('Comparison','fontsize',24);
        
        ylim([0.5, numel(meanData)+0.5]);
        yrange = ylim;
        
        
        switch factor_name
            case {'CPP'}
                xrange = [-0.3, 0.3];
            case {'RU'}
                xrange = [-0.6, 0.6];
            case {'Reward'}
                xrange = [-0.1, 0.1];
            case {'Residual'}
                xrange = [-0.3, 0.3];
        end
        xlim(xrange);
        
        
        for i = 1:numel(pvalData)
            
            idx_p = 0;
            if pvalData(i)<.05 & pvalData(i)>=0.01
                idx_p = 1;
                list_n = 0;
            elseif pvalData(i)<0.01 & pvalData(i)>=0.001
                idx_p = 1;
                list_n = [-0.5,0.5];
            elseif pvalData(i)<0.001
                idx_p = 1;
                list_n = [-1,0,1];
            end
            
            if idx_p==1
                y = i;
                if meanData(i)>=0
                    x = meanData(i)+semData(i) + 0.05*(max(xrange)-min(xrange));
                else
                    x = meanData(i)-semData(i) - 0.05*(max(xrange)-min(xrange));
                end
                
                list_n = list_n*0.25;
                for n = 1:numel(list_n)
                    newy = y + list_n(n) - 0.15;
                    text(x,newy,'*','HorizontalAlignment','left','VerticalAlignment', 'middle', 'fontsize',24);
                end
                
            end
            
        end
        
        title(factor_name, 'fontsize', 24);
        
        outputfile = fullfile(dirFig, sprintf('coef_effect_comparison_%s_subgraph%.2d_%s',coef_name, idx_subgraph, factor_name));
        print(outputfile,'-depsc','-opengl');
        
    end
end



% individual difference of normative learning and dynamic modulation
clear lesion_result;
list_factor = {'Dynamic modulation'};
nFactor = numel(list_factor);
normative_LR = allData.all.behav.behav_beta(:,3) + allData.all.behav.behav_beta(:,4);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(1, nComparison);
    lesion_result.(coef_name).sem = NaN(1, nComparison);
    lesion_result.(coef_name).pval = NaN(1, nComparison);
    
    % orignal data
    current_data.all = allData.all.reg.coef_behav.cpp_ru_vals.rel(:,1,:) +...
        allData.all.reg.coef_behav.cpp_ru_vals.rel(:,2,:);
    current_data.within = allData.within.reg.coef_behav.cpp_ru_vals.rel(:,1,:) +...
        allData.within.reg.coef_behav.cpp_ru_vals.rel(:,2,:);
    current_data.between = allData.between.reg.coef_behav.cpp_ru_vals.rel(:,1,:) +...
        allData.between.reg.coef_behav.cpp_ru_vals.rel(:,2,:);
    
    for idx_comparison = 1:nComparison
        
        x1 = current_data.(list_comparison{idx_comparison,1})(:,:,idx_subgraph);
        x2 = current_data.(list_comparison{idx_comparison,2})(:,:,idx_subgraph);
        
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
        
        lesion_result.(coef_name).mean(1,idx_comparison) = diff_mean;
        lesion_result.(coef_name).sem(1,idx_comparison) = diff_sem;
        lesion_result.(coef_name).pval(1,idx_comparison) = pval;
        
    end
    
end

% graph
plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
plotParameter.color.bar.edgecolor = {'none'};
plotParameter.color.error = [0, 0, 0];
plotParameter.linewidth = 2;

plotParameter.xticklabel = {'Within-All', 'Between-All', 'Within-Between'};

threshold_pval = .05;
for c = 1:nCoef
    
    coef_name = list_coef{c};
    
    current_mean = lesion_result.(coef_name).mean(1,:);
    current_sem = lesion_result.(coef_name).sem(1,:);
    current_pval = lesion_result.(coef_name).pval(1,:);
    
    %%%%% horizontal bar %%%%%
    for f = 1:nFactor
        
        factor_name = list_factor{f};
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*0.85;
        fg.position(4) = fg.position(4)*0.45;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        
        meanData = current_mean(f,:);
        semData = current_sem(f,:);
        pvalData = current_pval(f,:);
        
        hold on
        
        plot([0;0], [0.5;numel(meanData)+0.5],...
            'linestyle', '--',...
            'linewidth', 2,...
            'color', [0.5,0.5,0.5]);
        e = errorbar(meanData', [1:nComparison]', semData', 'horizontal');
        set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
        plot(meanData', [1:nComparison]',...
            'marker', 'o',...
            'linestyle', 'none',...
            'linewidth', 2,...
            'markersize', 12,...
            'markerfacecolor', [0.5,0.5,0.5],...
            'markeredgecolor', 'k');
        
        hold off
        
        set(gca,'YTick',[1:nComparison],'YTickLabel',list_type);
        set(gca, 'fontsize',18, 'linewidth',2);
        
        xlabel('Change of coefficient','fontsize',24);
        ylabel('Comparison','fontsize',24);
        
        ylim([0.5, numel(meanData)+0.5]);
        yrange = ylim;
        
        xrange = [-0.55, 0.55];
        xlim(xrange);
        
        
        for i = 1:numel(pvalData)
            
            idx_p = 0;
            if pvalData(i)<.05 & pvalData(i)>=0.01
                idx_p = 1;
                list_n = 0;
            elseif pvalData(i)<0.01 & pvalData(i)>=0.001
                idx_p = 1;
                list_n = [-0.5,0.5];
            elseif pvalData(i)<0.001
                idx_p = 1;
                list_n = [-1,0,1];
            end
            
            if idx_p==1
                y = i;
                if meanData(i)>=0
                    x = meanData(i)+semData(i) + 0.05*(max(xrange)-min(xrange));
                else
                    x = meanData(i)-semData(i) - 0.05*(max(xrange)-min(xrange));
                end
                
                list_n = list_n*0.25;
                for n = 1:numel(list_n)
                    newy = y + list_n(n) - 0.15;
                    text(x,newy,'*','HorizontalAlignment','left','VerticalAlignment', 'middle', 'fontsize',24);
                end
                
            end
            
        end
        
        title(factor_name, 'fontsize', 24);
        
        outputfile = fullfile(dirFig, sprintf('coef_corr_dynamic_comparison_%s_subgraph%.2d',coef_name, idx_subgraph));
        print(outputfile,'-depsc','-opengl');
        
    end
    
end


% individual difference of normative learning and average expression
clear lesion_result;
list_factor = {'Average expression'};
nFactor = numel(list_factor);
normative_LR = allData.all.behav.behav_beta(:,3) + allData.all.behav.behav_beta(:,4);
for c = 1:nCoef
    
    % coef_name
    coef_name = list_coef{c};
    lesion_result.(coef_name).mean = NaN(1, nComparison);
    lesion_result.(coef_name).sem = NaN(1, nComparison);
    lesion_result.(coef_name).pval = NaN(1, nComparison);
    
    % orignal data
    current_data.all = allData.all.subgraph.subject.(coef_name).mean;
    current_data.within = allData.within.subgraph.subject.(coef_name).mean;
    current_data.between = allData.between.subgraph.subject.(coef_name).mean;
    
    for idx_comparison = 1:nComparison
        
        x1 = current_data.(list_comparison{idx_comparison,1})(:,idx_subgraph);
        x2 = current_data.(list_comparison{idx_comparison,2})(:,idx_subgraph);
        
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
        
        lesion_result.(coef_name).mean(1,idx_comparison) = diff_mean;
        lesion_result.(coef_name).sem(1,idx_comparison) = diff_sem;
        lesion_result.(coef_name).pval(1,idx_comparison) = pval;
        
    end
    
end

% graph
plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
plotParameter.color.bar.edgecolor = {'none'};
plotParameter.color.error = [0, 0, 0];
plotParameter.linewidth = 2;

plotParameter.xticklabel = {'Within-All', 'Between-All', 'Within-Between'};

threshold_pval = .05;
for c = 1:nCoef
    
    coef_name = list_coef{c};
    
    current_mean = lesion_result.(coef_name).mean(1,:);
    current_sem = lesion_result.(coef_name).sem(1,:);
    current_pval = lesion_result.(coef_name).pval(1,:);
    
    %%%%% horizontal bar %%%%%
    for f = 1:nFactor
        
        factor_name = list_factor{f};
        
        figure;
        fg = fig_setting_default;
        fg.position = get(gcf, 'Position');
        fg.position(3) = fg.position(3)*0.85;
        fg.position(4) = fg.position(4)*0.45;
        set(gcf, 'Position', fg.position);
        set(gcf, 'PaperPosition', fg.position);
        
        
        meanData = current_mean(f,:);
        semData = current_sem(f,:);
        pvalData = current_pval(f,:);
        
        hold on
        
        plot([0;0], [0.5;numel(meanData)+0.5],...
            'linestyle', '--',...
            'linewidth', 2,...
            'color', [0.5,0.5,0.5]);
        e = errorbar(meanData', [1:nComparison]', semData', 'horizontal');
        set(e, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
        plot(meanData', [1:nComparison]',...
            'marker', 'o',...
            'linestyle', 'none',...
            'linewidth', 2,...
            'markersize', 12,...
            'markerfacecolor', [0.5,0.5,0.5],...
            'markeredgecolor', 'k');
        
        hold off
        
        set(gca,'YTick',[1:nComparison],'YTickLabel',list_type);
        set(gca, 'fontsize',18, 'linewidth',2);
        
        xlabel('Change of coefficient','fontsize',24);
        ylabel('Comparison','fontsize',24);
        
        ylim([0.5, numel(meanData)+0.5]);
        yrange = ylim;
        
        xrange = [-0.25, 0.25];
        xlim(xrange);
        
        
        for i = 1:numel(pvalData)
            
            idx_p = 0;
            if pvalData(i)<.05 & pvalData(i)>=0.01
                idx_p = 1;
                list_n = 0;
            elseif pvalData(i)<0.01 & pvalData(i)>=0.001
                idx_p = 1;
                list_n = [-0.5,0.5];
            elseif pvalData(i)<0.001
                idx_p = 1;
                list_n = [-1,0,1];
            end
            
            if idx_p==1
                y = i;
                if meanData(i)>=0
                    x = meanData(i)+semData(i) + 0.05*(max(xrange)-min(xrange));
                else
                    x = meanData(i)-semData(i) - 0.05*(max(xrange)-min(xrange));
                end
                
                list_n = list_n*0.25;
                for n = 1:numel(list_n)
                    newy = y + list_n(n) - 0.15;
                    text(x,newy,'*','HorizontalAlignment','left','VerticalAlignment', 'middle', 'fontsize',24);
                end
                
            end
            
        end
        
        title(factor_name, 'fontsize', 24);
        
        outputfile = fullfile(dirFig, sprintf('coef_corr_average_comparison_%s_subgraph%.2d',coef_name, idx_subgraph));
        print(outputfile,'-depsc','-opengl');
        
    end
    
end

