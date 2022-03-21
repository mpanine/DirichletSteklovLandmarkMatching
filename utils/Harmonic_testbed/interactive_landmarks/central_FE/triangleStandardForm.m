function [l1, tmax, a, b, Rloc, ind1, ind2] = triangleStandardForm(tri_inds, edge_lengths, landmark, desired_edge, R )
%TRIANGLESTANDARDFORM translates a triangle into standard form in order to
%compute the modified stiffness matrix
%tri_inds: indices of vertices of the desired triangle ex: [1 2 3]
%edge_lengths: lengths of the edges ( [1 2], [2 3] and [3 1] ) for eventual compatibility with intrinsic triangulations
%landmark: index of the currently studied landmark (the code does not support landmarks that are two edges or less apart)
%desired_edge: edge on which we are computing the stiffness
%R: internal radius of landmark (aka epsilon)

%OUTPUTS:
%l1, tmax, a, b: Usual paramters of standard triangle
%Rloc: radius at the zeroth vertex (R if landmark is part of the triangle, 0 otherwise)
%ind1, ind2: standard form indices of the vertices of desired_edge in (between 0 and 2)




if sum( tri_inds == landmark) == 1  %the landmark is part of the triangle

    while tri_inds(1) ~= landmark % position the landmark as the first vertex
        
        tri_inds = circshift( tri_inds, 1);
        edge_lengths = circshift(edge_lengths, 1); 
        
    end

    l1 = edge_lengths(1); %horizontal
    tmax = acos( ( edge_lengths(1)^2 + edge_lengths(3)^2 - edge_lengths(2)^2) / ( 2*edge_lengths(1)*edge_lengths(3) ) ); % law of cosines
    a = edge_lengths(3) * cos(tmax);
    b = edge_lengths(3) * sin(tmax);
    
    Rloc = R;
    
    ind1 = find(tri_inds == desired_edge(1) ) - 1; % -1 as need to re-index
    ind2 = find(tri_inds == desired_edge(2) ) - 1;
    
    

elseif sum( tri_inds == landmark) == 0 %the landmark is not part of the triangle (any order is fine)

    l1 = edge_lengths(1); %horizontal
    tmax = acos( ( edge_lengths(1)^2 + edge_lengths(3)^2 - edge_lengths(2)^2) / ( 2*edge_lengths(1)*edge_lengths(3) ) ); % law of cosines
    a = edge_lengths(3) * cos(tmax);
    b = edge_lengths(3) * sin(tmax);
    
    ind1 = find(tri_inds == desired_edge(1) ) - 1; % -1 as need to re-index
    ind2 = find(tri_inds == desired_edge(2) ) - 1;

    
    Rloc = 0.00001; %No internal radius.  TEMPORARY VERSION
    
else
    
   error('The triangle contains the same landmark a nonstandard number of times!') 
    
end
   


%% NORMALIZE THE LENGTHS

% 
% a = a/l1;
% b = b/l1;
% Rloc = Rloc/l1;
% l1 = 1;
% tmax - atan(b/a)
% tmax


% edge_lengths
% l1
% a
% b
% tmax

%%


if isempty(ind1)
    error('ind1 is empty: the desired edge vertex is not a summit of the triangle')
end

if isempty(ind2)
    error('ind2 is empty: the desired edge vertex is not a summit of the triangle')
end







end

