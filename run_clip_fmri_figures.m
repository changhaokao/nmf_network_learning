function run_clip_fmri_figures()

window_width = 10;
window_overlap = 8;


% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirFig = fullfile(dirType, fullfile(sprintf('figures_%d_%d', window_width, window_overlap),'brain'));
dirClip = fullfile(dirType, fullfile(sprintf('figures_%d_%d', window_width, window_overlap),'brain_clip'));
mkdir(dirClip);

% setting
temp_frame_x = [350:1200];
temp_frame_y = [160:1250];

type_align = 'square';

% list
list_map = {
    'NEURON_CPP_pos';
    'NEURON_CPP_neg';
    'NEURON_RU_pos';
    'NEURON_RU_neg';
    'NEURON_Rwd_pos';
    'NEURON_Rwd_neg';
    'subgraph04_pos';
    'subgraph04_neg';
    };

nMap = numel(list_map);

list_view = {
    'leftlateral';
    'leftmedial';
    'rightlateral';
    'rightmedial';
    };
nView = numel(list_view);

for m = 1:nMap
    
    % map
    map_name = list_map{m};
    
    clear matrix_A matrix_transparency;
    for v = 1:nView
        
        % view
        current_view = list_view{v};
        
        % load map
        inputfile = fullfile(dirFig, sprintf('%s_%s.png', map_name, current_view));
        outputfile = fullfile(dirClip, sprintf('clip_%s.png', map_name));
        
        [A,map,transparency] = imread(inputfile);
        transparency = transparency(temp_frame_x, temp_frame_y);
        A = A(temp_frame_x, temp_frame_y, :);
        
        [idx_x,idx_y] = find(transparency>0);
        frame_x = [min(idx_x):max(idx_x)];
        frame_y = [min(idx_y):max(idx_y)];
        transparency = transparency(frame_x, frame_y);
        A = A(frame_x, frame_y, :);
        
        % store
        matrix_A.(current_view) = A;
        matrix_transparency.(current_view) = transparency;
        
    end
    
    switch type_align
        case {'square'}
            A = [...
                [matrix_A.leftlateral, matrix_A.rightlateral];
                [matrix_A.leftmedial, matrix_A.rightmedial]];
            transparency = [...
                [matrix_transparency.leftlateral, matrix_transparency.rightlateral];
                [matrix_transparency.leftmedial, matrix_transparency.rightmedial]];
        case {'horizontal'}
            A = [matrix_A.leftlateral, matrix_A.leftmedial, matrix_A.rightmedial, matrix_A.rightlateral];
            transparency = [matrix_transparency.leftlateral, matrix_transparency.leftmedial, matrix_transparency.rightmedial, matrix_transparency.rightlateral];
    end
    
    imwrite(A, outputfile, 'alpha', transparency);
    
end

