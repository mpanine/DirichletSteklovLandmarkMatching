function eF = extendSegmentFunctions(F, segment, nv)
% pads the functions contained in F to nv vertices using nans

% F: matrix of functions on segment. size: length(segment) * numfun * extra_dim (Used in Dirichlet-Steklov)
% segment: list of vertices in the segment
% nv: number of desired vertices

% 

if length(size(F))==2

    numf = size(F,2);
    eF = nan(nv, numf);

    eF(segment,:) = F;
    
else %There is an additional dimension used in the Dirichlet-Steklov basis.
    
    eF = nan(nv, size(F,2), size(F,3));    
    eF(segment,:,:) = F;
    
end

end