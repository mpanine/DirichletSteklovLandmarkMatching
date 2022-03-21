function [S] = compute(S_in, shape_settings, cache_settings)
%COMPUTE Computes a Shape object
%   Input: S_in
%       - a mesh structure with surface (X,Y,Z, TRIV)
%       - the mesh filename
%       - shape_settings: the settings object for the shape
%       - cache_settings: the settings object for the cache

if ischar(S_in)
    S = shape.io.read_shape(S_in);
    filepath = S_in(1:end-length(S.name));
elseif isstruct(S_in)
    S = S_in;
    if ~isfield(S,'nv'), S.nv = length(S.surface.X); end
    if ~isfield(S,'nf'), S.nf = size(S.surface.TRIV,1); end
    filepath ='/';
else
    error('Unsupported input type.')
end

S.settings = shape_settings;

cached_shape = cache.load(cache_settings, S.name, S.settings);

if length(fieldnames(cached_shape)) == 0
%     [S.W, S.A, S.area, S.sqrt_area] = shape.geometry.compute_area_weights(S);

%     [S.geometry.E2V, S.geometry.T2E, S.geometry.E2T, S.geometry.T2T] = ...
%     shape.geometry.connectivity(S.surface.TRIV);

    S.geometry.normal = cross(...
    S.surface.VERT(S.surface.TRIV(:,1),:) - S.surface.VERT(S.surface.TRIV(:,2),:),...
    S.surface.VERT(S.surface.TRIV(:,1),:) - S.surface.VERT(S.surface.TRIV(:,3),:)...
                              );
    S.area = sqrt(sum(S.geometry.normal.^2, 2))/2;
    S.geometry.normal = S.geometry.normal./repmat(sqrt(sum(S.geometry.normal.^2, 2)), [1, 3]);
    S.sqrt_area = sqrt(sum(S.area));

%     S.geometry.SqEdgeLength = sum(...
% (S.surface.VERT(S.geometry.E2V(:,1),:) - S.surface.VERT(S.geometry.E2V(:,2),:)).^2 ...
%                                   ,2);

%     [S.W,~] = shape.geometry.per_edge_laplacian(abs(S.geometry.T2E),...
%                                  S.geometry.E2V,...
%                                  S.area, S.geometry.SqEdgeLength);
%     S.W = -S.W;
%     [S.A,~] = shape.geometry.per_edge_area(S.geometry.E2T, S.geometry.E2V, S.area);


    [S.W, S.A, S.area, S.sqrt_area] = shape.geometry.compute_area_weights( S );
    S.A = S.A/trace(S.A);
    
    S.D = sparse(diag(sum(S.A,1)));
%     S.W = -S.W;
    
    % compute the Laplacian basis
    if S.settings.computeMHB
        S.laplacianBasis = basis.compute_MHB(S,shape_settings.MHB_settings);
    end
    
    % get the one-ring neighbor
    if S.settings.findNeigh
        S.vtx_neigh = shape.geometry.find_one_ring_neigh(S);
    end
    
    % edge list and weight
    if S.settings.findEdge
        [S.Elist,S.EdgeWeight] = shape.geometry.get_edge_list(S);
    end
    
    % compute the geodesic distance
    if S.settings.computeGeoDist
        if exist([cache_settings.directory,'Gammas/',S.name,'.mat'],'file')
            Gamma_file = load([cache_settings.directory,'Gammas/',S.name,'.mat']);
            Gamma = Gamma_file.Gamma;
        else
            if ~exist([cache_settings.directory,'Gammas/'],'dir')
                mkdir([cache_settings.directory,'Gammas/']);
            end
            Gamma = shape.geometry.compute_geodesic_dist_matrix(S);
            save([cache_settings.directory,'Gammas/',S.name,'.mat'], 'Gamma');
        end
        S.Gamma = Gamma;
    end
    
    % compute the .obj file
    if S.settings.computeObj
        objpath = [cache_settings.directory,'Objs/',S.name,'.obj'];
        if ~exist(objpath,'file')
            if ~exist([cache_settings.directory,'Objs/'],'dir')
                mkdir([cache_settings.directory,'Objs/']);
            end
            shape.io.writeObj(objpath, S);
        end
        S.objpath = objpath;
    end
    
    % Normals and areas
    if S.settings.computeNormals
        [S.normals_vtx, S.normals_face] = shape.geometry.compute_vtx_and_face_normals(S);
    end
    cache.save(S, S.name, cache_settings);
else
    S = cached_shape.object;
end
end
