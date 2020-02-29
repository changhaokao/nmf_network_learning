function run_surfice_cmd_subgraph()

window_width = 10;
window_overlap = 8;

% directory
dirType = 'result_24mc_CSF_WM/result_raw';
dirSurfice = '/Applications/Surfice';
dirMap = fullfile(dirType, sprintf('result_subsystem_%d_%d', window_width, window_overlap));
dirFig = fullfile(dirType, sprintf('figures_%d_%d', window_width, window_overlap),'brain');
mkdir(dirFig);

% app
exe = fullfile(dirSurfice, 'surfice.app/Contents/MacOS/surfice');

% mesh
list_hemisphere = {'left', 'right'};
nHemisphere = numel(list_hemisphere);
mesh.left = 'BrainMesh_ICBM152Left.mz3';
mesh.right = 'BrainMesh_ICBM152Right.mz3';

% view
list_view = {'lateral', 'medial'};
nView = numel(list_view);

% NEURON
list_map = {
        {'node_strength_subgraph04_pos.nii.gz', ''}, 'subgraph04_pos';
        {'', 'node_strength_subgraph04_neg.nii.gz'}, 'subgraph04_neg';
    };


nMap = size(list_map,1);
list_color = {'Red-Yellow', 'Blue-Green'};
zmin = 0.5;
zmax = 1;


% surfice
for m = 1:nMap
    
    current_map = list_map{m,1};
    nEffect = numel(current_map);
    
    for h = 1:nHemisphere
        
        current_hemisphere = list_hemisphere{h};
        current_mesh = mesh.(current_hemisphere);
        switch current_hemisphere
            case {'left'}
                list_angle = [270, 90]; % lateral, medial
            case {'right'}
                list_angle = [90, 270]; % lateral, medial
        end
        
        cmd = sprintf('%s -S', exe);
        cmd = sprintf('%s "begin', cmd);
        cmd = sprintf('%s meshload(''%s'');', cmd, current_mesh);
        
        idx_effect = 0;
        for e = 1:nEffect
            
            filename = current_map{e};
            if isempty(filename)
                continue
            end
            inputfile = fullfile(dirMap, filename);
            data = load_nii(inputfile);
            if ~any(data.img(:)>0)
                continue
            end
            idx_effect = idx_effect + 1;
            cmd = sprintf('%s overlayload(''%s'');', cmd, inputfile);
            cmd = sprintf('%s overlaycolorname(%d, ''%s'');', cmd, idx_effect, list_color{e});
            cmd = sprintf('%s overlayminmax(%d, %.1f, %.1f);', cmd, idx_effect, zmin, zmax);
            
        end
        
        for v = 1:nView
            
            current_view = list_view{v};
            current_angle = list_angle(v);
            cmd = sprintf('%s azimuthelevation(%d,0);', cmd, current_angle);
            
            filename = sprintf('%s_%s%s', list_map{m,2}, current_hemisphere, current_view);
            outputfile = fullfile(dirFig, filename);
            cmd = sprintf('%s savebmp(''%s'');', cmd, outputfile);
            
        end % end of view
        
        cmd = sprintf('%s quit;', cmd);
        cmd = sprintf('%s end."', cmd);
        system(cmd);
        
    end % end of hemisphere
end % end of map


