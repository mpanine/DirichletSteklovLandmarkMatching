function [ to_keep ] = remove_border_landmarks_slow( SHAPE, landmarks )
    %REMOVE_BORDER_LANDMARKS Summary of this function goes here
    %   Detailed explanation goes here
    
    %% V0 = TOO COSTLY (neighbor computation)
    W = SHAPE.W;
    % Exclude the vertices in this region from the set of acceptable landmarks
    bd = unique(reshape(calc_boundary_edges(SHAPE.surface.TRIV),[],1));
    exclude_neighs = bd;
    for idx=bd
        % 2-ring neighs of indices
        neighs = find(W(idx,:));
        exclude_neighs = [exclude_neighs,find(ismember(landmarks,neighs))];
        neighs2 = find(sum(W(neighs,:),1));
        exclude_neighs = [exclude_neighs,find(ismember(landmarks,neighs2))];
    end
    exclude_neighs = unique(exclude_neighs);

    to_keep = find(~ismember(landmarks,exclude_neighs));
end

