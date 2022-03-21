function [ shapePairList ] = generate_random_shapePairList( datasetDir, dataset, length )
    %GENERATE_RANDOM_SHAPEPAIRLIST Generates a random set of shape pairs of
    %length "length" for a dataset.
    dataDir = [datasetDir,dataset,'/'];
    
    switch dataset
        case {'Faust', 'Shrec07_Fourleg', 'Shrec16_Partial/holes', 'Shrec16_Partial/cuts'}
            files = dir([dataDir,'*.off']);
            assert(size(files, 1)>2, 'Selected directory is empty.');
        case 'Shrec19'
            files = dir([dataDir,'*.obj']);
            assert(size(files, 1)>2, 'Selected directory is empty.');
        case {'Faust_remeshed', 'SCAPE_remeshed', 'TOSCA_nonIsometric'}
            files = dir([dataDir,'vtx_5k/*.off']);
            assert(size(files, 1)>2, 'Selected directory is empty.');
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
    end
    
    switch dataset
        case 'Shrec16_Partial/holes'
            files_holes = dir([dataDir,'*.off']);
            start_file_idx_holes = shape.datasetTools.find_start_file_idx(files_holes);
            tar_fidx = [];
            for fi=start_file_idx_holes:size(files_holes, 1)
                fname = [files_holes(fi).name];
                [~, name, ~] = fileparts(fname);
                if isempty(strfind(name,'_null'))
                    tar_fidx(end+1) = fi;
                end
            end
            
            assert(size(tar_fidx,2) > length, 'length is greater than the numher of available files.');
            
            rand_order = randperm(size(tar_fidx,2), length);
            tar_fidx = tar_fidx(rand_order);
            shapePairList = {};
            for i=1:length
                fname = [files_holes(tar_fidx(i)).name];
                [~, shapeName, ~] = fileparts(fname);
                src_name = [shapeName,'_null'];
                tar_name = shapeName;

                shapePairList{end+1} = [src_name, ' ', tar_name];
            end
        case 'Shrec16_Partial/cuts'
            files_cuts = dir([dataDir,'*.off']);
            start_file_idx_cuts = shape.datasetTools.find_start_file_idx(files_cuts);
            tar_fidx = [];
            for fi=start_file_idx_cuts:size(files_cuts, 1)
                fname = [files_cuts(fi).name];
                [~, name, ~] = fileparts(fname);
                if isempty(strfind(name,'_null'))
                    tar_fidx(end+1) = fi;
                end
            end
            
            assert(size(tar_fidx,2) > length, 'length is greater than the numher of available files.');
            
            rand_order = randperm(size(tar_fidx,2), length);
            tar_fidx = tar_fidx(rand_order);
            shapePairList = {};
            for i=1:length
                fname = [files_cuts(tar_fidx(i)).name];
                [~, shapeName, ~] = fileparts(fname);
                src_name = [shapeName,'_null'];
                tar_name = shapeName;

                shapePairList{end+1} = [src_name, ' ', tar_name];
            end
        case 'Shrec19'
            pair_list = load(sprintf('%sPAIRS_list_SHREC19_connectivity.txt', dataDir));
            pair_list = pair_list(1:208,:); % We don't have landmarks after 26
            rand_indices = randperm(size(pair_list, 1), length);
            shapePairList = {};
            for i=1:size(rand_indices,2)
                src_name = pair_list(rand_indices(i), 1);
                tar_name = pair_list(rand_indices(i), 2);
                
                shapePairList{end+1} = [char(string(src_name)), ' ', char(string(tar_name))];
            end
        case 'Shrec19_remeshed'
            pair_list = load(sprintf('%svtx_5k/remeshed_rescaled/PAIRS_list_SHREC19_connectivity.txt', dataDir));
            rand_indices = randperm(size(pair_list, 1), length);
            
            shapePairList = {};
            for i=1:size(rand_indices,2)
                src_name = num2str(pair_list(rand_indices(i), 1));
                tar_name = num2str(pair_list(rand_indices(i), 2));
                
                shapePairList{end+1} = [src_name, ' ', tar_name];
            end
        case 'TOSCA_Isometric'
            shapePairList = {};
            for i=1:length
                shape_type = randi([1, size(files_list, 2)],1, 1);
                files = files_list{shape_type};
                shapes = randperm(size(files, 1),2);
                src_name = regexprep(files(shapes(1)).name,'.off','');
                tar_name = regexprep(files(shapes(2)).name,'.off','');
                shapePairList{end+1} = [src_name, ' ', tar_name];
            end
        otherwise
            m = size(files, 1);
            k = randperm(m/2*(m-1),length);
            q = max(floor(sqrt(8*(k-1) + 1)/2 + 3/2), 2);
            p = max(k - (q-1).*(q-2)/2, 2);
    
            selected_indices = sortrows([p;q]',[2 1]);

            shapePairList = {};
            for i=1:size(selected_indices,1)
                src_name = regexprep(files(selected_indices(i,1)).name,'.off','');
                tar_name = regexprep(files(selected_indices(i,2)).name,'.off','');

                shapePairList{end+1} = [src_name, ' ', tar_name];
            end
    end
end

