function group_rsquared_within_between_effect()

window_width = 10;
window_overlap = 8;

% directory
dirType = fullfile('result_24mc_CSF_WM','result_raw');
dirFig = fullfile(dirType, sprintf('figures_%d_%d', window_width, window_overlap));
mkdir(dirFig);


% setting
nPerm = 10000;
nSubgraph = 10;


%%%%% within-subject: task effects %%%%%
allReg = group_behav_subgraph_coef('', 0);
rsquared.within_subject = mean(allReg.reg.rsquared.cpp_ru_vals.rel)';

%%%%% between-subject %%%%%
y = allReg.behav.behav_beta(:,3)+allReg.behav.behav_beta(:,4);
x1 = allReg.subgraph.subject.rel.mean;
x2 = allReg.reg.coef_behav.cpp_ru_vals.rel(:,1,:)+allReg.reg.coef_behav.cpp_ru_vals.rel(:,2,:);
x2 = reshape(x2,32,[]);

rsquared.between_subject = NaN(nSubgraph, 1);
for i = 1:nSubgraph
    
    x = [x1(:,i),x2(:,i)];
    [beta, dev, stats] = glmfit(x, y);
    rsquared.between_subject(i,1) = 1-var(stats.resid)./var(y);
    
end




%%%%%%%%%%%%%%%%%
%%%%% graph %%%%%
%%%%%%%%%%%%%%%%%
list_pair = {
    'within_subject', 'between_subject';
    };
nPair = size(list_pair,1);

for i = 1:nPair
    
    name_pair1 = list_pair{i,1};
    name_pair2 = list_pair{i,2};
    
    figure;
    fig_setting_default;
    
    x = rsquared.(name_pair1);
    y = rsquared.(name_pair2);
    
    [rho,p] = permutation_corr(x,y,nPerm);
    
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
    
    
    xlabel(sprintf('Within-subject effects'),'fontsize',24);
    xrange = [0.01, 0.018];
    
    ylabel(sprintf('Between-subject effects'),'fontsize',24);
    yrange = [-0.1, 0.35];
    
    xlim(xrange);
    ylim(yrange);
    
    xpos = xrange(1) + (xrange(2)-xrange(1))*0.05;
    ypos = yrange(1) + (yrange(2)-yrange(1))*0.9;
    
    if p>=.001
        text_p = sprintf('p=%.3f',p);
    elseif p<.001 & p>=.0001
        text_p = 'p<.001';
    elseif p<.0001
        text_p = 'p<.0001';
    end
    
    text(xpos,ypos,sprintf('r=%.3f\n%s', rho, text_p), 'fontsize',24);
    
    outputfile = fullfile(dirFig, sprintf('corr_subgraph_%s_%s', name_pair1, name_pair2));
    print(outputfile,'-depsc','-opengl');
    
end


