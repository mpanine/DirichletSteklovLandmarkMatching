function [ output_args ] = generate_landmark_files( dataset, datasetDir, saveDir )
    %GENERATE_LANDMARK_FILES Generates consistent landmarks for all
    %shapes within a given dataset, stored to disk.
    
    shape_settings = shape.Settings();
    cache_settings = cache.Settings(dataset);
    
    dataDir = [datasetDir,dataset,'/',];
    
    if nargin < 3
        saveDir = datasetDir;
    end
    
    switch dataset
        case 'Shrec16_Partial'
            null_models_dir = [dataDir,'/null/'];
            holes_models_dir = [dataDir,'/holes/'];
            cuts_models_dir = [dataDir,'/cuts/'];
            
            files_null = dir([null_models_dir,'*.off']);
            start_file_idx_null = shape.datasetTools.find_start_file_idx(files_null);
            files_holes = dir([holes_models_dir,'*.off']);
            start_file_idx_holes = shape.datasetTools.find_start_file_idx(files_holes);
            files_cuts = dir([cuts_models_dir,'*.off']);
            start_file_idx_cuts = shape.datasetTools.find_start_file_idx(files_cuts);
        case 'Shrec19'
            files = dir([dataDir, '*.obj']);
        case 'Shrec19_remeshed'
            files = dir([dataDir, 'vtx_5k/remeshed_rescaled/*.obj']);
        case {'Faust_remeshed','SCAPE_remeshed'}
            files = dir([dataDir, 'vtx_5k/*.off']);
            start_file_idx = 1;
            first_mesh_name = files(start_file_idx).name;
        case 'TOSCA_nonIsometric'
            files = dir([dataDir, 'vtx_5k/*.off']);
            start_file_idx = 1;
            first_mesh_name = files(start_file_idx).name;
        case 'TOSCA_Isometric'
            files_cat = dir([dataDir, 'vtx_5k/cat*.off']);
            files_centaur = dir([dataDir, 'vtx_5k/centaur*.off']);
            files_david = dir([dataDir, 'vtx_5k/david*.off']);
            files_dog = dir([dataDir, 'vtx_5k/dog*.off']);
            files_horse = dir([dataDir, 'vtx_5k/horse*.off']);
            files_michael = dir([dataDir, 'vtx_5k/michael*.off']);
            files_victoria = dir([dataDir, 'vtx_5k/victoria*.off']);
            files_wolf = dir([dataDir, 'vtx_5k/wolf*.off']);
            files_list = {files_cat,...
                     files_centaur,...
                     files_david,...
                     files_dog,...
                     files_horse,...
                     files_michael,...
                     files_victoria,...
                     files_wolf};
        otherwise
            files = dir([dataDir, '*.off']);
            start_file_idx = shape.datasetTools.find_start_file_idx(files);
            first_mesh_name = files(start_file_idx).name;
    end
    
    if ~exist([saveDir,dataset,'/Landmarks/'], 'dir')
       mkdir([saveDir,dataset,'/Landmarks/']);
    end
    
    switch dataset
        case {'Faust_690', 'Faust'}
            first_shape = shape.compute([dataDir,first_mesh_name], shape_settings, cache_settings);
            initial_samples = shape.geometry.euclidean_fps(first_shape, 20);
            remaining_samples = 1:first_shape.nv;
            remaining_samples(initial_samples) = [];
            
            for i=start_file_idx:size(files, 1)
                mesh_name = regexprep(files(i).name,'.off','');
                
                % concatenating fps-computed landmarks at 20 locations and
                % all remaining points
                landmarks = [initial_samples', remaining_samples];
                save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
                clear landmarks;
            end
        case 'Shrec07_Fourleg'
%             full_ordered_landmarks = [1,5,10,13,15, 18, 2,9,6,7,11, 19, 3,8, 17, 4,12, 21, 14,16, 20, 21];
            
            for i=start_file_idx:size(files, 1)
                if ~isempty(strfind(files(i).name, '.off'))
                    mesh_name = regexprep(files(i).name,'.off','');

                    % Loading the .vts file associated with the mesh
                    vts = load([dataDir,mesh_name,'.vts']);
                    % Retrieving the landmarks (first column of the .vts
                    % array)
                    landmarks = vts(:, 1)+1;
%                     landmarks = landmarks(full_ordered_landmarks);
                    save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
                    clear landmarks;
                end
            end
        case {'Faust_remeshed','SCAPE_remeshed'}
            first_shape = shape.compute([dataDir,'vtx_5k/',first_mesh_name], shape_settings, cache_settings);
            mesh_name = regexprep(first_mesh_name,'.off','');
            
            full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
            
            landmarks_Src = shape.geometry.euclidean_fps(first_shape, first_shape.nv,1);
            [~,permutation] = ismember(landmarks_Src,full_landmarks);
            permutation = permutation(find(permutation));
            landmarks = full_landmarks(permutation);
            save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
            clear landmarks;
            clear full_landmarks;
            
            
            for i=start_file_idx+1:size(files, 1)
                mesh_name = regexprep(files(i).name,'.off','');
                
                full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
                landmarks = full_landmarks(permutation);
                save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
                clear landmarks;
                clear full_landmarks;
            end
%             first_shape = shape.compute([dataDir,'vtx_5k/',first_mesh_name], shape_settings, cache_settings);
%             mesh_name = regexprep(first_mesh_name,'.off','');
%             
%             full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
%             sampling_shape.surface.VERT = first_shape.surface.VERT(full_landmarks, :);
%             
%             initial_samples = shape.geometry.euclidean_fps(sampling_shape, 20);
%             remaining_samples = full_landmarks;
%             remaining_samples(initial_samples) = [];
%             % concatenating fps-computed landmarks at 20 locations and
%                 % all remaining points
%             landmarks_Src = [initial_samples; remaining_samples]';
%             landmarks = full_landmarks(landmarks_Src);
%             save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
%             clear landmarks;
%             clear full_landmarks;
%             
%             
%             for i=start_file_idx+1:size(files, 1)
%                 mesh_name = regexprep(files(i).name,'.off','');
%                 
%                 full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
%                 landmarks = full_landmarks(landmarks_Src);
%                 save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
%                 clear landmarks;
%                 clear full_landmarks;
%             end
        case 'TOSCA_Isometric'
            for f_idx=1:size(files_list, 2)
                files = files_list{f_idx};
                first_mesh_name = files(1).name;
                start_file_idx = 1;
                first_shape = shape.compute([dataDir,'vtx_5k/',first_mesh_name], shape_settings, cache_settings);
                mesh_name = regexprep(first_mesh_name,'.off','');

                full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
                sampling_shape.surface.VERT = first_shape.surface.VERT(full_landmarks, :);

                initial_samples = shape.geometry.euclidean_fps(sampling_shape, 20);
                remaining_samples = full_landmarks;
                remaining_samples(initial_samples) = [];
                % concatenating fps-computed landmarks at 20 locations and
                    % all remaining points
                landmarks_Src = [initial_samples; remaining_samples]';
                landmarks = full_landmarks(landmarks_Src);
                save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
                clear landmarks;
                clear full_landmarks;

                for i=start_file_idx+1:size(files, 1)
                    mesh_name = regexprep(files(i).name,'.off','');

                    full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
                    landmarks = full_landmarks(landmarks_Src);
                    save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
                    clear landmarks;
                    clear full_landmarks;
                end
            end
            
        case 'TOSCA_nonIsometric'
            first_shape = shape.compute([dataDir,'vtx_5k/',first_mesh_name], shape_settings, cache_settings);
            mesh_name = regexprep(first_mesh_name,'.off','');
            
            sampling_vert_idx = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
            sampling_shape.surface.VERT = first_shape.surface.VERT(sampling_vert_idx, :);
            
            initial_samples = shape.geometry.euclidean_fps(sampling_shape, 20);
            remaining_samples = sampling_vert_idx;
            remaining_samples(initial_samples) = [];
            % concatenating fps-computed landmarks at 20 locations and
                % all remaining points
            landmarks = [initial_samples; remaining_samples];
            
            
            
            save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
            landmarks_ref = knn(sampling_vert_idx, landmarks);%knnsearch(sampling_vert_idx, landmarks);
            clear landmarks;
            
            
            for i=start_file_idx+1:size(files, 1)
                mesh_name = regexprep(files(i).name,'.off','');
                
                full_landmarks = load([dataDir,'/vtx_5k/corres/',mesh_name,'.vts']);
                landmarks = full_landmarks(landmarks_ref);
                save([saveDir,dataset,'/Landmarks/',mesh_name,'.mat'],'landmarks');
                clear landmarks;
                clear full_landmarks;
            end
            
        case 'Shrec19'
            pair_list = load('Data/Shrec19/PAIRS_list_SHREC19_connectivity.txt');
            matches_dir = [ dataDir, 'matches/FARMgt_txt/' ];
            
            for i=1:size(files, 1)
                i
                Src = shape.compute(sprintf('%s/%d.obj', dataDir,i), shape_settings, cache_settings);
                
                initial_samples = shape.geometry.euclidean_fps(Src, 20);
                remaining_samples = 1:Src.nv;
                remaining_samples(initial_samples) = [];
                landmarks = [initial_samples', remaining_samples];
                
                save(sprintf('%s%s/Landmarks/%d.mat', saveDir, dataset, i),'landmarks');
                landmarks_Src = landmarks;
                clear landmarks;
                
                Tar_list = pair_list(pair_list(:, 1) == i, 2);
                
                for ls_idx=1:size(Tar_list, 1)
                    j = Tar_list(ls_idx)
                    pair_corres = load(sprintf('%s%d_%d.txt', matches_dir, i, j));
                    landmarks = pair_corres(landmarks_Src);
                    save(sprintf('%s%s/Landmarks/%d_%d.mat', saveDir, dataset, i, j),'landmarks');
                    clear landmarks;
                end
            end
            
        case 'Shrec19_remeshed'
            pair_list = load('Data/Shrec19_remeshed/vtx_5k/remeshed_rescaled/PAIRS_list_SHREC19_connectivity.txt');
            matches_dir = [ dataDir, 'vtx_5k/remeshed_rescaled/groundtruth/' ];
            
            for i=1:size(files, 1)
                i
                Src = shape.compute(sprintf('%svtx_5k/remeshed_rescaled/%d.obj', dataDir,i), shape_settings, cache_settings);
                
                initial_samples = shape.geometry.euclidean_fps(Src, 20);
                remaining_samples = 1:Src.nv;
                remaining_samples(initial_samples) = [];
                landmarks = [initial_samples', remaining_samples];
                
                save(sprintf('%s%s/Landmarks/%d.mat', saveDir, dataset, i),'landmarks');
                landmarks_Src = landmarks;
                clear landmarks;
                
                Tar_list = pair_list(pair_list(:, 1) == i, 2);
                
                for ls_idx=1:size(Tar_list, 1)
                    j = Tar_list(ls_idx)
                    pair_corres = load(sprintf('%s%d_%d.map', matches_dir, i, j));
                    landmarks = pair_corres(landmarks_Src);
                    save(sprintf('%s%s/Landmarks/%d_%d.mat', saveDir, dataset, i, j),'landmarks');
                    clear landmarks;
                end
            end
        case 'Shrec16_Partial'
            % Holes
            
            if ~exist([saveDir,dataset,'/Landmarks/holes/'], 'dir')
               mkdir([saveDir,dataset,'/Landmarks/holes/']);
            end
            for fi=1:size(files_holes, 1)
                fname = [files_holes(fi).name];
                [~, name, ~] = fileparts(fname);
                
                cat_parts = strsplit(name,'_');
                
                if strcmp(cat_parts{2},'cat')&&strcmp(cat_parts{4},'-1')
                    if isempty(strfind(name,'_null'))
                        null_model = regexp(name,'[^_0-9]*','match');
                        null_model = [null_model{2} ''];

                        null_shape = shape.compute([null_models_dir, null_model], shape_settings, cache_settings);
                        holes_shape = shape.compute([holes_models_dir, name], shape_settings, cache_settings);

                        data = load([holes_models_dir, name,'.baryc_gt']);
                        % Removing correspondences that do not exist (for cat.off
                        % shapes, due to the non-manifold faces)
                        data(find(data(:,1)>size(null_shape.surface.TRIV,1)),:) = [];
                        all_GT_correspondences = bar2idx(data, null_shape.surface.TRIV);

                        clear landmarks;
                        % Generating landmarks on the partial model via FPS
                        % sampling
                        landmarks = shape.geometry.euclidean_fps(holes_shape, 200, 1);
                        save([saveDir,dataset,'/Landmarks/holes/', name,'.mat'], 'landmarks');
                        % Retrieving the corresponding vertices on the null shape
                        landmarks = all_GT_correspondences(landmarks);
                        save([saveDir,dataset,'/Landmarks/holes/', name,'_null.mat'], 'landmarks');
                        clear landmarks;
                    end
                end
            end
            
            % Cuts
            
            if ~exist([saveDir,dataset,'/Landmarks/cuts/'], 'dir')
               mkdir([saveDir,dataset,'/Landmarks/cuts/']);
            end
            for fi=1:size(files_cuts, 1)
                fname = [files_cuts(fi).name];
                [~, name, ~] = fileparts(fname);
                
                cat_parts = strsplit(name,'_');
                
                if strcmp(cat_parts{2},'cat')&&strcmp(cat_parts{4},'1')
                    if isempty(strfind(name,'_null'))
                        null_model = regexp(name,'[^_0-9]*','match');
                        null_model = [null_model{2} ''];

                        null_shape = shape.compute([null_models_dir, null_model], shape_settings, cache_settings);
                        cuts_shape = shape.compute([cuts_models_dir, name], shape_settings, cache_settings);

                        data = load([cuts_models_dir, name,'.baryc_gt']);
                        % Removing correspondences that do not exist (for cat.off
                        % shapes, due to the non-manifold faces)
                        data(find(data(:,1)>size(null_shape.surface.TRIV,1)),:) = [];
                        all_GT_correspondences = bar2idx(data, null_shape.surface.TRIV);

                        clear landmarks;
                        % Generating landmarks on the partial model via FPS
                        % sampling
                        landmarks = shape.geometry.euclidean_fps(cuts_shape, 200, 1);
                        save([saveDir,dataset,'/Landmarks/cuts/', name,'.mat'], 'landmarks');
                        % Retrieving the corresponding vertices on the null shape
                        max_idx = size(all_GT_correspondences,1);
                        if max(landmarks)>max_idx
                                good_lmks = landmarks(landmarks<=max_idx);
                                good_lmks = all_GT_correspondences(good_lmks);
                                bad_lmks = landmarks(landmarks>max_idx);
                                remaining = setdiff(all_GT_correspondences,good_lmks);
                                bad_lmks = remaining(1:size(bad_lmks,1));
                                landmarks = [good_lmks; bad_lmks];
                        else
                            landmarks = all_GT_correspondences(landmarks);
                        end
                        save([saveDir,dataset,'/Landmarks/cuts/', name,'_null.mat'], 'landmarks');
                        clear landmarks;
                    end
                end
            end
        otherwise
            throw(MException('GENERATELANDMARKFILES:datasetUnknown',...
                'The dataset is not handled. Please add it to the switch in the constructor of shape.datasetTools.generate_landmark_files.'));
    end
    
end
