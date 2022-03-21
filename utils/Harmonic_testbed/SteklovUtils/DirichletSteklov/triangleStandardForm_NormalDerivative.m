function tri_out = triangleStandardForm_NormalDerivative(tri, edge)
%TRIANGLESTANDARDFORM_REF puts a triangle in standard form with respect to a
%given edge for normal derivative calculation purposes (The edge of interst is between vertices 1 and 2)
%tri: 3 indices defining a triangle
%edge: 2 indices defining an edge
%If the edge is not contained in the triangle, the triangle in unchanged.


tri_out = tri;
fledge = flip(edge);


for i = 0:2
    
    t = circshift(tri_out, i); %Current form of the triangle.
    
    if sum(t(1:2) == edge) == 2
        
        tri_out = t; % The desired edge is in the correct position
        break;
        
    elseif sum(t(1:2) == fledge) == 2
        
        tri_out = t; % The desired flipped edge is in the correct position
        break;
        
    end
    
end





end

