function [idx] = geodesic_fps(dist_matrix,k,seed)
%GEODESIC_FPS Summary of this function goes here
%   Detailed explanation goes here
    % uniqueness of vertex list (required by distances method)
    nv= size(dist_matrix,1);
    
    if(nargin<4)
        idx = randi(nv,1);
    else
        idx = seed;
    end
    dists = dist_matrix(:,idx);
    for i = 1:k
        maxi = find(dists == max(dists));
        maxi = maxi(1);
        idx = [idx; maxi];
        newdists = dist_matrix(:,maxi);
        dists = min(dists,newdists);
    end

    idx = idx(2:end);
end

