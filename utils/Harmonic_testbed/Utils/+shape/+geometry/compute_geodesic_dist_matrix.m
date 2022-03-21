function D = compute_geodesic_dist_matrix(S, samples, useLegacy)
% Compute the Geodesic distance matrix of a mesh (or a set of samples)
    M = S.surface;
    M.n = S.nv;

    if nargin==1 || isempty(samples)
        samples = 1:M.n;
    end
    
    if nargin<3
        useLegacy = false;
    end

    fprintf('Compute geodesic distance matrix for %d samples...',length(samples)); tic;
    if useLegacy
        % Legacy code using fastmarchmex    
        march = fastmarchmex('init', int32(M.TRIV-1), double(M.VERT(:,1)), double(M.VERT(:,2)), double(M.VERT(:,3)));

        D = zeros(length(samples));

        for i=1:length(samples)
        %     fprintf('(%d/%d)\n', i, length(samples));
            source = inf(M.n,1);
            source(samples(i)) = 0;
            d = fastmarchmex('march', march, double(source));
            D(:,i) = d(samples);
        end

        fastmarchmex('deinit', march);

        D = 0.5*(D+D');
    else
        X2 = sum(M.VERT.^2,2);
        distance = sqrt(repmat(X2,1,M.n)+repmat(X2',M.n,1)-2*M.VERT*M.VERT');

        if ~isfield(S,'Elist')
            [S.Elist,S.EdgeWeight] = shape.geometry.get_edge_list(S);
        end
        adj_mat = zeros(M.n);
        adj_mat(sub2ind(size(adj_mat),S.Elist(:,1),S.Elist(:,2))) = 1;
        G = sparse(adj_mat.*distance);
        G = G+G';

        D = dijkstra(G, samples);
    end
    t = toc; fprintf('done:%.4fs\n',t);
end
