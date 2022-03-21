function [err,curve,err_raw,curve_raw] = compute_errors_and_curves_SHREC19(M,thresh,c_matches,c_gt_matches)
%COMPUTE_RECON_ERRORS_AND_CURVES Summary of this function goes here
%   Detailed explanation goes here
% uniqueness of vertex list (required by distances method)
[src_idx, src_ia, src_ic]=unique(c_matches);
[sink_idx, sink_ia, sink_ic]=unique(c_gt_matches);

% Compute dijkstra over graph
geo_d = distances(M.G,src_idx,sink_idx);

% come back to 6890 x 6890 of distance between FARM_gt and method
dist_matrix = geo_d(src_ic,sink_ic);

% if matches are on disconnected components, we use euclidean
% metric

[ii, jj] = find(isinf(dist_matrix));
if not(isempty(ii))
    eucl=vecnorm(M.surface.VERT(c_matches(ii),:)-M.surface.VERT(c_gt_matches(jj),:),2,2);
    idx = sub2ind(size(dist_matrix),ii,jj);
    dist_matrix(idx)=eucl;
end

% Normalize error
err = diag(dist_matrix);
diam = max(max(dist_matrix));

err_raw = err;
err = err ./ diam; 

% Compute Curve
curve_raw = calc_err_curve(err_raw, [0:0.001:0.1]);

curve = calc_err_curve(err, thresh);
end

