function ND = normalDerivativeAtLandmarkBoundary(surface, F, boundary, boundary_mass, landmark, num_unrefined_triangles)
%normalDerivativeAtBoundary: Computes the normal derivative of F at the
%boundary. ASSUMES CIRCULAR BOUNDARIES.

%surface: structure containing fields VERT and TRIV
%boundary: list of boundary vertices in the standard order (neighbors are part of edges)
%boundary_mass: Steklov mass matrix for the relevant boundary

num_vert = length(boundary); %number of vertices and also of triangles in a circular boundary

%% Finding the relevant triangles.
 

RelevantTRIV = zeros(num_vert, 3); % Triangles used in the computation. Will be put in standard order.

for t = 1:num_vert % Finding the relevant triangles and putting them in standard order.
    
    edge = [boundary(t) boundary( mod(t,num_vert) + 1 ) ];    
    tri_indices = findTrianglesContainingEdge_restricted(surface.TRIV, edge, num_unrefined_triangles+1, size(surface.TRIV,1)); %findTrianglesContainingEdge(surface, edge);
    
    for i = 1:length(tri_indices)
        
        if sum(surface.TRIV(tri_indices(i),:) == landmark) == 0 %The landmark is not in the triangle
           
            current_tri = surface.TRIV(tri_indices(i),:);
            
        end
                
    end
    
    RelevantTRIV(t,:) = triangleStandardForm_NormalDerivative(current_tri, edge);  
    
end


%% Finding the relevant geometric quantities

L1 = sqrt( sum( ( surface.VERT(RelevantTRIV(:,1),:) - surface.VERT(RelevantTRIV(:,2),:) ).^2 , 2) ); %Edge lengths
L2 = sqrt( sum( ( surface.VERT(RelevantTRIV(:,2),:) - surface.VERT(RelevantTRIV(:,3),:) ).^2 , 2) );
L3 = sqrt( sum( ( surface.VERT(RelevantTRIV(:,3),:) - surface.VERT(RelevantTRIV(:,1),:) ).^2 , 2) );

Cosines = (L3.^2 + L1.^2 - L2.^2) ./ (2*L1.*L3);
Angles = acos(Cosines);
Sines = sin(Angles);

Cx = L3 .* Cosines; % "x" and "y" positions of the non-boundary point in the standard form.
Cy = L3 .* Sines;
 
%% Computing the normal derivatives at the edges

ND_edges = ( F(RelevantTRIV(:,3)) - F(RelevantTRIV(:,1)) - F(RelevantTRIV(:,2)) .* Cx ./ L1 ) ./ Cy;


%% Computing the normal derivatives at the vertices

ND = ( ND_edges .* L1  +  circshift( ND_edges .* L1, 1) ) ./ ( L1 + circshift(L1, 1) );


%% Normalizing according to the mass matrix

nn = sqrt( ND' * boundary_mass * ND );

ND = ND / nn;



end