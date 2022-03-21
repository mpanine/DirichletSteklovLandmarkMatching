function [vert_idx] = bar2idx(data, TRIV)
% Convert barycentric matches to vertex index matches
    [~,tri_vert_idx] = max(data(:,2:end),[],2);
    
    vert_idx = TRIV(sub2ind(size(TRIV),data(:,1),tri_vert_idx));
end