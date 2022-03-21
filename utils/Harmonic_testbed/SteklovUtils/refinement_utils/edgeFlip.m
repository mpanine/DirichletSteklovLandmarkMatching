function surface_out = edgeFlip(surface, edge)
%EDGEFLIP flips the prescribed edge and outputs the updated surface
%surface: structure with fields TRIV, X, Y, Z and VERT
%edge: edge to be flipped

tri_indices = findTrianglesContainingEdge(surface, edge); %finds the relevant triangles



if length(tri_indices) == 2 % There is an edge to flip 
    
    %TODO: FIX POTENTIAL TRIANGLE FLIPS AND OTHER SUCH ISSUES
    
    surface_out = surface;
    
    % Put the triangles in standard form: the edge of interest comes last (between vertices 3 and 1) -----------------------
    tri_1 = triangleStandardForm_Ref( surface.TRIV(tri_indices(1),:), edge);
    tri_2 = triangleStandardForm_Ref( surface.TRIV(tri_indices(2),:), edge);
    
    surface_out.TRIV(tri_indices(1), :) = [tri_2(2) tri_2(3) tri_1(2)]; % Flip the edge
    surface_out.TRIV(tri_indices(2), :) = [tri_1(2) tri_1(3) tri_2(2)];
    
%     tri_1(3) == tri_2(1)
%     tri_2(3) == tri_1(1)
    
else % The edge does not correspond to exactly two triangles: do nothing
    
    surface_out = surface;   
    fprintf(sprintf('Found %d edge adjacent triangles instead of the expected 2\n',length(tri_indices)))
    
end




end

