function [errs] = calc_distortion_err(matches,gt_matches, D)
%CALC_GEO_ERR Summary of this function goes here
%   Detailed explanation goes here
    nm = length(matches);
    errs = zeros(nm,1);
    
    for i=1:nm
        errs(i) = D(matches(i),gt_matches(i));
    end
    
    errs = errs ./ max(max(D));
end

