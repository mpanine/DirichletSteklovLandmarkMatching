function eF = extendSegmentFunctionsIntegrals(F, F_bound, Mass_Boundary, segment, landmarks, nv)
% Extends the functions contained in F (undefined at the landmarks) to the landmarks.
% The values at the landmarks are given by boundary integrals.


% F: matrix of functions on segment. size: length(segment) * numfun * extra_dim (Used in Dirichlet-Steklov)
% F_bound: F restricted to the boundary elements
% Mass_Boundary: Mass boundary 
% segment: list of vertices in the segment (refined mesh without the landmarks)
% nv: number of desired vertices

% 

if length(size(F))==2

    numf = size(F,2);
    eF = zeros(nv, numf);

    eF(segment,:) = F;
    
else %There is an additional dimension used in the Dirichlet-Steklov basis.
    
    eF = zeros(nv, size(F,2), size(F,3));    
    eF(segment,:,:) = F;
    
end


if length(size(F))==2
    
    for i = 1:length(landmarks)
    
        eF(landmarks(i),:) = sum(Mass_Boundary{i}*F_bound{i},1);% / sum(sum( Mass_Boundary{i} )); % Values at the landmarks are integrals over the boundary.
    
    end
        
elseif length(size(F))==3
    
    for i = 1:length(landmarks) %Cycle over the boundaries with nonzero integrals
        
        eF(landmarks(i),:,i) = sum(Mass_Boundary{i}*F_bound{i},1);% / sum(sum( Mass_Boundary{i} )); % Values at the landmarks are integrals over the boundary.

    end
    
end
    
    













end