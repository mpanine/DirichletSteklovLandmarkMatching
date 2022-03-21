function [ to_keep ] = remove_border_landmarks( SHAPE, landmarks )
    %REMOVE_BORDER_LANDMARKS Summary of this function goes here
    %   Detailed explanation goes here
    
%     %% V0 = TOO COSTLY (neighbor computation)
%     W = SHAPE.W;
%     % compute the 1-ring neighbor
%     get_neighbor = @(W) cellfun(@(x,i) setdiff(reshape(find(x),[],1),i),...
%         mat2cell(W,ones(size(W,1),1),size(W,2)),num2cell(1:size(W,1))','UniformOutput',false);
%     vtx_neigh = get_neighbor(W);
%     % Exclude the vertices in this region from the set of acceptable landmarks
%     bd = unique(reshape(calc_boundary_edges(SHAPE.surface.TRIV),[],1));
%     exclude_neighs = bd;
%     for idx=exclude_neighs
%         size(exclude_neighs)
%         neigh = vertcat(vtx_neigh{idx});
%         exclude_neighs = [exclude_neighs; neigh];
%     end
%     exclude_neighs = unique(exclude_neighs);


    edges = SHAPE.Elist;
    vertices = SHAPE.surface.VERT;
    triangles = SHAPE.surface.TRIV;
    
    % find border vertices
    exclude_neighs = unique(reshape(calc_boundary_edges(triangles),[],1));
    
    % find max edge length along border
    c1 = find(ismember(edges(:,1),exclude_neighs));
    c2 = find(ismember(edges(:,2),exclude_neighs));
    c = [c1; c2];
    edge_lengths = sqrt(sum((vertices(edges(c,1),:)-vertices(edges(c,2),:)).^2,2));
    max_elen = max(edge_lengths);

    if max_elen>0
        [indices,~] = nearestneighbour(exclude_neighs',vertices', 'Radius', max_elen);
        exclude_neighs = unique(indices(:));
        exclude_neighs(exclude_neighs==0) = [];

        % indices to keep among landmarks
        % --> filter correct landmarks by applying:
        % landmarks = landmarks(to_keep);
    end
    to_keep = find(~ismember(landmarks,exclude_neighs));
end

