function allReg = group_behav_subgraph_coef(typename, idx_fig, result_type, window_setting, idx_interaction_remove)

if nargin<1
    typename = '';
end
if nargin<2
    idx_fig = 1;
end

if nargin<3
    result_type = 'result_raw';
    
    % result_type = 'result_glm_bold_temporal_derivative';
    
end

if nargin<4
    window_setting.window_width = 10;
    window_setting.window_overlap = 8;
end
window_width = window_setting.window_width;
window_overlap = window_setting.window_overlap;



% directory
dirType = fullfile('result_24mc_CSF_WM',result_type);

dirSubgraph = fullfile(dirType, sprintf('result_consensus%s_%d_%d', typename, window_width, window_overlap));
switch typename
    case {'_remove_interaction_individual'}
        dirSubgraph = fullfile(dirSubgraph, sprintf('remove_%d', idx_interaction_remove));
        idx_fig = 0;
end

dirSubsystem = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));
dirHazard = 'behavior_mle';

if idx_fig==1
    dirFig = fullfile(dirType, sprintf('figures%s_%d_%d', typename, window_width, window_overlap));
    mkdir(dirFig);
end


% load behav_param
load('behav_param.mat');
load('behav_model_slope.mat');

% load subgraph
load(fullfile(dirSubgraph,'result_ens.mat'));

% setting
nSubgraph = size(coef_ens,1);

nRun = 4;
nTRperRun = 226;
TR_unused_initial = [1:6];
TR_unused_end = [nTRperRun-(window_width-2):nTRperRun];

nTR_unused = numel(TR_unused_initial) + numel(TR_unused_end);
TR = 2.5; % second
nTRperRun_include = nTRperRun-nTR_unused;

nWindow = (nTRperRun-numel(TR_unused_initial)-window_width)/(window_width-window_overlap) + 1;

nbin = 20;

nPerm = 10000;

list_coef = {'rel'};

% sublist
sublist = textread('sublist','%s');
nSub = numel(sublist);

% sort subject data
cpp = [];
ru = [];
vals = [];
edge = [];
noise = [];
update = [];
LR = [];
PE = [];
isCP = [];
coins = [];

currentPrediction = [];
yHat = [];

outcomeTimes = [];

subj_index = [];
block_index = [];
blockNum = [];
badtrial = [];


behav_beta = NaN(size(upSlope));
for s = 1:nSub
    
    subname = sublist{s};
    idx_sub = strcmp(behav_param.subjName,subname);
    
    idx_include = (idx_sub);
    
    subjNum = behav_param.subjNum(idx_include);
    subjNum = unique(subjNum);
    
    behav_beta(s,:) = upSlope(subjNum,:);
    
    outcomeTimes = [
        outcomeTimes;
        behav_param.outcomeTimes(idx_include)
        ];
    cpp = [
        cpp;
        behav_param.modpCha(idx_include)
        ];
    ru = [
        ru;
        behav_param.modRelUnc(idx_include)
        ];
    vals = [
        vals;
        behav_param.vals(idx_include)
        ];
    edge = [
        edge;
        behav_param.edge(idx_include)
        ];
    noise = [
        noise;
        behav_param.noise(idx_include)
        ];
    update = [
        update;
        behav_param.update(idx_include)
        ];
    LR = [
        LR;
        behav_param.LR(idx_include)
        ];
    PE = [
        PE;
        behav_param.PE(idx_include)
        ];
    isCP = [
        isCP;
        behav_param.isChangeTrial(idx_include)
        ];
    coins = [
        coins;
        behav_param.coinsCought(idx_include)
        ];
    subj_index = [
        subj_index;
        ones(sum(idx_include),1)*s
        ];
    block_index = [
        block_index;
        behav_param.newBlock(idx_include)
        ];
    blockNum = [
        blockNum;
        behav_param.blockNum(idx_include)
        ];
    badtrial = [
        badtrial;
        behav_param.badtrial(idx_include)
        ];
    yHat = [
        yHat;
        behav_param.yHat(idx_include)
        ];
end
block_index = cumsum(block_index);

residual = [abs(update)-abs(yHat)];
residual = residual./abs(PE);
residual(residual>=1) = 1;
residual(residual<=-1) = -1;

idx_time = (outcomeTimes>(max(TR_unused_initial)*TR)) & (outcomeTimes<((min(TR_unused_end)-1)*TR));
idx_include = (badtrial~=1)&idx_time&(~isnan(LR));
outcomeTimes = outcomeTimes(idx_include);

cpp = cpp(idx_include);
ru = ru(idx_include);
vals = vals(idx_include);
noise = noise(idx_include);
update = update(idx_include);
LR = LR(idx_include);
PE = PE(idx_include);
residual = residual(idx_include);
residual_abs = abs(residual);

subj_index = subj_index(idx_include);
blockNum = blockNum(idx_include);


coef_ens = coef_ens';
nTime = size(coef_ens,1);
coef_pos = coef_ens(1:nTime/2,:);
coef_neg = coef_ens(nTime/2+1:nTime,:);
coef_pos_trial = [];
coef_neg_trial = [];
window_relative_rms_trial = [];
for s = 1:nSub
    
    idx_sub = (subj_index==s);
    
    for b = 1:nRun
        
        idx_block = (blockNum==b);
        idx_select = idx_sub&idx_block;
        
        current_outcomeTimes = outcomeTimes(idx_select);
        current_outcomeTimes = current_outcomeTimes - max(TR_unused_initial)*TR;
        behav_TR_index = current_outcomeTimes/TR + 1;
        
        coef_TR_index = (window_width-window_overlap)*([1:nWindow]'-1) + 1;
        idx_period = [1:nWindow] + ((b-1)+(s-1)*nRun)*nWindow;
        
        coef_pos_block = coef_pos(idx_period,:);
        coef_neg_block = coef_neg(idx_period,:);
        coef_pos_block = interp1(coef_TR_index,coef_pos_block,behav_TR_index);
        coef_neg_block = interp1(coef_TR_index,coef_neg_block,behav_TR_index);
        
        coef_pos_trial = [coef_pos_trial; coef_pos_block];
        coef_neg_trial = [coef_neg_trial; coef_neg_block];
        
    end
    
end

coef_pos = coef_pos_trial;
coef_neg = coef_neg_trial;
coef_rel = coef_pos - coef_neg;

idx_include = ~isnan(sum(coef_rel,2));
cpp = cpp(idx_include);
ru = ru(idx_include);
vals = vals(idx_include);
noise = noise(idx_include);
update = update(idx_include);
LR = LR(idx_include);
PE = PE(idx_include);
residual = residual(idx_include);
residual_abs = abs(residual);

subj_index = subj_index(idx_include);
blockNum = blockNum(idx_include);

coef_pos = coef_pos(idx_include,:);
coef_neg = coef_neg(idx_include,:);
coef_rel = coef_rel(idx_include,:);

allReg.behav.behav_beta = behav_beta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% regression: coef~behav_var %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpp, ru, vals
nCoef = numel(list_coef);

design_matrix = [cpp,ru,vals,residual];

nReg = size(design_matrix,2);

for c = 1:nCoef
    
    name_coef = list_coef{c};
    
    allReg.reg.coef_behav.cpp_ru_vals.(name_coef) = NaN(nSub,nReg,nSubgraph);
    
    eval(sprintf('subgraph_coef = coef_%s;',name_coef));
    
    for i = 1:nSubgraph
        
        for s = 1:nSub
            
            idx_sub = (subj_index==s);
            
            x = design_matrix(idx_sub,:);
            y = subgraph_coef(idx_sub,i);
            
            x = zscore(x);
            
            [beta, dev, stats] = glmfit(x, y);
            
            r_squared = 1-var(stats.resid)./var(y);
            
            allReg.reg.coef_behav.cpp_ru_vals.(name_coef)(s,:,i) = beta(2:end);
            allReg.reg.rsquared.cpp_ru_vals.(name_coef)(s,i) = r_squared;
            
        end
        
        
    end
    
end




%%%%% graph %%%%%
if idx_fig
    plotParameter.color.bar.facecolor = {[0.5,0.5,0.5]};
    plotParameter.color.bar.edgecolor = {'none'};
    plotParameter.color.error = [0, 0, 0];
    plotParameter.linewidth = 2;
    plotParameter.xticklabel = {'CPP', 'RU', 'Reward', 'Residual'};
    plotParameter.ylabel = 'Coefficient';
    plotParameter.barwidth = 0.5;
    
    for c = 1:nCoef
        
        name_coef = list_coef{c};
        
        for s = 1:nSubgraph
            
            figure;
            fg = fig_setting_default();
            
            subjData = allReg.reg.coef_behav.cpp_ru_vals.(name_coef)(:,:,s);
            nSub = size(subjData,1);
            
            groupMean = mean(subjData,1);
            groupSEM = std(subjData,0,1)./sqrt(nSub);
            
            meanData = groupMean;
            semData = groupSEM;
            
            % bar graph
            hold on
            
            plot([0;numel(meanData)+1], [0;0],...
                'Color', [0.5,0.5,0.5],...
                'LineWidth', 2,...
                'LineStyle', '-');
            
            h = bar(meanData);
            barsx = h.XData;
            e = errorbar(barsx, meanData, semData,...
                'Marker', 'none',...
                'LineStyle', 'none'...
                );
            set(h,...
                'FaceColor', plotParameter.color.bar.facecolor{1},...
                'EdgeColor', plotParameter.color.bar.edgecolor{1},...
                'barwidth', plotParameter.barwidth,...
                'LineWidth', plotParameter.linewidth...
                );
            set(e,...
                'Color', plotParameter.color.error,...
                'LineWidth', plotParameter.linewidth...
                );
            hold off
            
            hold on
            xMargin = 0.04;
            yMargin = 0.04;
            for i = 1:size(subjData,2)
                
                current_data = subjData(:,i);
                [xOffset] = smartJitter(current_data,xMargin,yMargin);
                xpos = i+0.35+xOffset;
                plot(xpos, current_data,...
                    'linestyle', 'none',...
                    'marker', 'o',...
                    'markersize', 12,...
                    'markerfacecolor', 'w',...
                    'markeredgecolor', 'k',...
                    'linewidth', 2);
                
            end
            
            hold off
            
            set(gca,...
                'XTick', barsx,...
                'XTickLabel', plotParameter.xticklabel...
                );
            
            set(gca,'fontsize',24,'linewidth',2);
            
            ylabel(plotParameter.ylabel,...
                'fontsize', 32);
            
            xlim([0.5, numel(groupMean)+0.5]);
            yrange = ylim;
            yrange = yrange + [-0.2,0.2].*(yrange(2)-yrange(1));
            ylim(yrange);
            
            
            [h,p,ci, stats] = ttest(subjData);
            for i = 1:numel(p)
                
                idx_p = 0;
                if p(i)<.05 & p(i)>=0.01
                    idx_p = 1;
                    text_p = '*';
                elseif p(i)<0.01 & p(i)>=0.001
                    idx_p = 1;
                    text_p = '**';
                elseif p(i)<0.001
                    idx_p = 1;
                    text_p = '***';
                end
                
                if idx_p==1
                    x = i;
                    if groupMean(i)>=0
                        y = groupMean(i)+groupSEM(i) + 0.05*(max(yrange)-min(yrange));
                    else
                        y = groupMean(i)-groupSEM(i) - 0.05*(max(yrange)-min(yrange));
                    end
                    text(x,y,text_p,'HorizontalAlignment','center','fontsize',32);
                    
                end
                
            end
            
            title(sprintf('Subgraph %d', s), 'fontsize',32);
            
            outputfile = fullfile(dirFig, sprintf('reg_coef_behav_cpp_ru_vals_%s_subgraph%.2d', name_coef,s));
            print(outputfile,'-depsc','-opengl');
            
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% property of subgraph %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all
within_connection = NaN(nSubgraph,1);
between_connection = NaN(nSubgraph,1);
for s = 1:nSubgraph
    
    % load file
    filename = fullfile(dirSubsystem, sprintf('subgraph_%.2d_subnet_subsystem.mat',s));
    load(filename);
    
    val_diag = diag(subsystem_subnet);
    within_connection(s) = mean(val_diag);
    
    val_off_diag = triu(subsystem_subnet,1);
    idx_include = logical(triu(ones(size(subsystem_subnet)),1));
    val_off_diag = val_off_diag(idx_include);
    between_connection(s) = mean(val_off_diag);
    
end
allReg.subgraph.all.within_connection = within_connection;
allReg.subgraph.all.between_connection = between_connection;
allReg.subgraph.all.within_ratio = (within_connection-between_connection)./(within_connection+between_connection);


nCoef = numel(list_coef);
subgraph_mean = NaN(nSubgraph,nCoef);
for c = 1:nCoef
    
    name_coef = list_coef{c};
    eval(sprintf('subgraph_coef = coef_%s;',name_coef));
    
    for s = 1:nSubgraph
        
        subgraph_mean(:,c) = mean(subgraph_coef,1);
        
    end
    
end
allReg.subgraph.all.mean = subgraph_mean;

% subject
nCoef = numel(list_coef);
for c = 1:nCoef
    
    name_coef = list_coef{c};
    
    eval(sprintf('subgraph_coef = coef_%s;',name_coef));
    
    allReg.subgraph.subject.(name_coef).mean = NaN(nSub,nSubgraph);
    
    for s = 1:nSub
        
        idx_sub = (subj_index==s);
        
        current_data = subgraph_coef(idx_sub,:);
        
        allReg.subgraph.subject.(name_coef).mean(s,:) = mean(current_data,1);
        
    end
    
end


%%%%% graph %%%%%
if idx_fig
    
    nCoef = numel(list_coef);
    
    list_pair = {
        'mean', 'within_ratio';
        'mean', 'within_connection';
        'mean', 'between_connection';
        };
    nPair = size(list_pair,1);
    
    for i = 1:nPair
        
        name_pair1 = list_pair{i,1};
        name_pair2 = list_pair{i,2};
        
        for c = 1:nCoef
            
            name_coef = list_coef{c};
            
            figure;
            fig_setting_default;
            
            y = allReg.subgraph.all.(name_pair1)(:,c);
            x = allReg.subgraph.all.(name_pair2);
            
            switch typename
                case {''}
                    [rho,p] = permutation_corr(x,y,nPerm);
                otherwise
                    [rho,p] = corr(x, y);
            end
            
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
                'linewidth', 2,...
                'Marker', 'o',...
                'MarkerSize', 24,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'w');
            
            hold off
            
            for j = 1:nSubgraph
                text(x(j), y(j), num2str(j), 'HorizontalAlignment', 'Center', 'fontsize', 18);
            end
            set(gca,'fontsize',24,'linewidth',2);
            
            ylabel('Average expression','fontsize',24);
            
            switch name_pair2
                case {'within_ratio'}
                    xlabel('Relative strength','fontsize',24);
                case {'within_connection'}
                    xlabel('Within-system strength','fontsize',24);
                case {'between_connection'}
                    xlabel('Between-system strength','fontsize',24);
            end
            
            xrange = xlim;
            yrange = ylim;
            switch name_pair2
                case {'within_ratio'}
                    xrange = [-0.7, 0.7];
                    yrange = [-10, 20];
                    xlim(xrange);
                    ylim(yrange);
                    xpos = xrange(1) + (xrange(2)-xrange(1))*0.05;
                    ypos = yrange(1) + (yrange(2)-yrange(1))*0.9;
                    set(gca,'XTick',[-0.6,-0.4,-0.2,0,0.2,0.4,0.6]);
                case {'within_connection'}
                    xrange = [0, 0.01];
                    yrange = [-10, 20];
                    xlim(xrange);
                    ylim(yrange);
                    xpos = xrange(1) + (xrange(2)-xrange(1))*0.05;
                    ypos = yrange(1) + (yrange(2)-yrange(1))*0.9;
                    set(gca,'XTick', [0:0.005:0.01], 'XTickLabel', [0:0.005:0.01]);
                case {'between_connection'}
                    xrange = [0.0025, 0.0055];
                    yrange = [-10, 20];
                    xlim(xrange);
                    ylim(yrange);
                    xpos = xrange(1) + (xrange(2)-xrange(1))*0.75;
                    ypos = yrange(1) + (yrange(2)-yrange(1))*0.9;
                    set(gca,'XTick', [0.0025:0.001:0.0055], 'XTickLabel', [0.0025:0.001:0.0055]);
            end
            
            if p>=.001
                text_p = sprintf('p=%.3f',p);
            elseif p<.001 & p>=.0001
                text_p = 'p<.001';
            elseif p<.0001
                text_p = 'p<.0001';
            end
            
            text(xpos,ypos,sprintf('r=%.3f\n%s', rho, text_p), 'fontsize',24);
            
            outputfile = fullfile(dirFig, sprintf('corr_subgraph_%s_%s_%s', name_pair1, name_coef, name_pair2));
            print(outputfile,'-depsc','-opengl');
            
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% average expression ~ behav_reg (CPP+RU) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normative_LR = allReg.behav.behav_beta(:,3) + allReg.behav.behav_beta(:,4);
if idx_fig
    
    switch result_type
        case {'result_raw'}
            location_ratio.rel = {
                [0.78, 0.9];
                [0.78, 0.9];
                [0.78, 0.15];
                [0.78, 0.15];
                [0.78, 0.15];
                [0.78, 0.9];
                [0.78, 0.9];
                [0.78, 0.9];
                [0.78, 0.15];
                [0.78, 0.9];
                };
        case {'result_glm_bold_temporal_derivative'}
            location_ratio.rel = {
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.15];
                [0.75, 0.15];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.15];
                [0.75, 0.9];
                };
    end
    
    
    list_measure = {'mean'};
    nMeasure = numel(list_measure);
    
    nCoef = numel(list_coef);
    
    for c = 1:nCoef
        
        name_coef = list_coef{c};
        
        for m = 1:nMeasure
            
            name_measure = list_measure{m};
            
            for s = 1:nSubgraph
                
                figure;
                fg = fig_setting_default();
                
                x = normative_LR;
                y = allReg.subgraph.subject.(name_coef).(name_measure)(:,s);
                
                switch typename
                    case {''}
                        [rho,p] = permutation_corr(x,y,nPerm);
                    otherwise
                        [rho,p] = corr(x, y);
                end
                
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
                    'linewidth', 2,...
                    'Marker', 'o',...
                    'MarkerSize', 12,...
                    'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', [0.5,0.5,0.5]);
                
                hold off
                
                
                set(gca,'fontsize',32,'linewidth',2);
                
                xlabel('Normative learning', 'fontsize', 32);
                
                ylabel('Average expression', 'fontsize', 32);
                title(sprintf('Subgraph %d',s), 'fontsize',32);
                
                xrange = [-0.7, 2.7];
                xlim(xrange);
                set(gca, 'XTick', [-0.5:0.5:2.5]);
                
                yrange = ylim;
                yrange = yrange + [-0.05,0.05].*(yrange(2)-yrange(1));
                ylim(yrange);
                
                xpos = xrange(1)+(xrange(2)-xrange(1))*location_ratio.rel{s}(1);
                ypos = yrange(1)+(yrange(2)-yrange(1))*location_ratio.rel{s}(2);
                
                if p>=.001
                    text_p = sprintf('p=%.3f',p);
                elseif p<.001 & p>=.0001
                    text_p = 'p<.001';
                elseif p<.0001
                    text_p = 'p<.0001';
                end
                
                text(xpos,ypos,sprintf('r=%.3f\n%s', rho, text_p), 'fontsize',24);
                
                outputfile = fullfile(dirFig, sprintf('corr_normative_LR_subgraph_%s_%s_subgraph%.2d', name_coef, name_measure,s));
                print(outputfile,'-depsc','-opengl');
                
            end
            
        end
        
    end
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% beta_subgraph ~ behav_reg (CPP+RU) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normative_LR = allReg.behav.behav_beta(:,3) + allReg.behav.behav_beta(:,4);
if idx_fig
    
    switch result_type
        case {'result_raw'}
            location_ratio.rel = {
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.15];
                [0.75, 0.15];
                [0.75, 0.15];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.15];
                };
        case {'result_glm_bold_temporal_derivative'}
            location_ratio.rel = {
                [0.75, 0.9];
                [0.75, 0.9];
                [0.05, 0.15];
                [0.75, 0.15];
                [0.05, 0.15];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                [0.75, 0.9];
                };
    end
    
    nCoef = numel(list_coef);
    
    for c = 1:nCoef
        
        name_coef = list_coef{c};
        
        for s = 1:nSubgraph
            
            beta_subgraph = allReg.reg.coef_behav.cpp_ru_vals.(name_coef)(:,:,s);
            
            figure;
            fg = fig_setting_default();
            
            x = normative_LR;
            y = beta_subgraph(:,1)+beta_subgraph(:,2);
            
            switch typename
                case {''}
                    [rho,p] = permutation_corr(x,y,nPerm);
                otherwise
                    [rho,p] = corr(x, y);
            end
            
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
                'linewidth', 2,...
                'Marker', 'o',...
                'MarkerSize', 12,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', [0.5,0.5,0.5]);
            
            hold off
            
            
            set(gca,'fontsize',32,'linewidth',2);
            
            xlabel('Normative learning', 'fontsize', 32);
            
            
            ylabel('Dynamic modulation', 'fontsize', 32);
            title(sprintf('Subgraph %d',s), 'fontsize',32);
            
            xrange = [-0.7, 2.7];
            xlim(xrange);
            set(gca, 'XTick', [-0.5:0.5:2.5]);
            
            yrange = ylim;
            
            xpos = xrange(1)+(xrange(2)-xrange(1))*location_ratio.rel{s}(1);
            ypos = yrange(1)+(yrange(2)-yrange(1))*location_ratio.rel{s}(2);
            
            if p>=.001
                text_p = sprintf('p=%.3f',p);
            elseif p<.001 & p>=.0001
                text_p = 'p<.001';
            elseif p<.0001
                text_p = 'p<.0001';
            end
            
            text(xpos,ypos,sprintf('r=%.3f\n%s', rho, text_p), 'fontsize',24);
            
            outputfile = fullfile(dirFig, sprintf('corr_normative_LR_subgraph_modulation_%s_subgraph%.2d', name_coef, s));
            print(outputfile,'-depsc','-opengl');
            
        end
        
    end
    
end



