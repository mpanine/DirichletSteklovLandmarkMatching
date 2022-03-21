function [  ] = Shrec16_Partial_generate_pair_shape( datasetDir )
    %SHREC16_PARTIAL_GENERATE_PAIR_SHAPE Generates redundant .off mesh
    %files for the null meshes in sthe Shrec16_Partial dataset, so that the
    %caching mechanism can store their configuration for the different 
    %holes and cuts models.
    
    dataDir = [datasetDir,'Shrec16_Partial','/',];
    
    null_models_dir = [dataDir,'/null/'];
    holes_models_dir = [dataDir,'/holes/'];
    cuts_models_dir = [dataDir,'/cuts/'];
    
    files_holes = dir([holes_models_dir,'*.off']);
    start_file_idx_holes = find_start_file_idx(files_holes);
    files_cuts = dir([cuts_models_dir,'*.off']);
    start_file_idx_cuts = find_start_file_idx(files_cuts);
    
    for fi=start_file_idx_holes:size(files_holes, 1)
        fname = [files_holes(fi).name];
        [~, name, ~] = fileparts(fname);
        
        fprintf(sprintf('Processing %s\n',name));

        null_model = regexp(name,'[^_0-9]*','match');
        null_model = [null_model{2} ''];
        
        S_null = shape.io.read_shape([null_models_dir, null_model, '.off']);
        
        shape.io.writeOFF([holes_models_dir, name, '_null.off'], S_null.surface.VERT, S_null.surface.TRIV);
    end
    
    for fi=start_file_idx_cuts:size(files_cuts, 1)
        fname = [files_cuts(fi).name];
        [~, name, ~] = fileparts(fname);
        
        fprintf(sprintf('Processing %s\n',name));

        null_model = regexp(name,'[^_0-9]*','match');
        null_model = [null_model{2} ''];
        
        S_null = shape.io.read_shape([null_models_dir, null_model, '.off']);
        
        shape.io.writeOFF([cuts_models_dir, name, '_null.off'], S_null.surface.VERT, S_null.surface.TRIV);
    end
end

