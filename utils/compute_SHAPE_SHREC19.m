function [ S ] = compute_SHAPE_SHREC19( filePath )
    %COMPUTE_MESH_NORMALIZED Summary computes the mesh with all elements required for
    %the wavelet computations (also for all baselines).
    % boundary_type = 'neumann' or 'dirichlet'
    % /!\ Although the vertex positions are divided by sqrt_area, the value
    % of sqrt_area remains the original value (i.e. is not set to 1) /!\
    
    shrec19_data = load([filePath '.mat']);
    N.surface.VERT = shrec19_data.Shape_df.VERT;
    N.surface.TRIV = shrec19_data.Shape_df.TRIV;
    N.G = shrec19_data.G;
    N.name = num2str(shrec19_data.name);%convertCharsToStrings(shrec19_data.name);
    
    
    
    N.nv = shrec19_data.Shape_df.n;
    N.nf = shrec19_data.Shape_df.m;
    
    N.surface.m = N.nf;
    N.surface.n = N.nv;
    N.surface.X = N.surface.VERT(:,1);
    N.surface.Y = N.surface.VERT(:,2);
    N.surface.Z = N.surface.VERT(:,3);
    
%     N = shape.compute(N,shape_settings,cache_settings);
    [N.W, N.A, N.area, N.sqrt_area] = shape.geometry.compute_area_weights( N );
    N.A = N.A/trace(N.A);
    
    N.D = sparse(diag(sum(N.A,1)));
    
    
    [N.Elist,~] = shape.geometry.get_edge_list(N);
    
    
    
    S = struct;
    S.landmarks = shrec19_data.smpl_matches;
    S.SHAPE = N;
end



