function tri_out = triangleStandardForm_Ref(tri, edge)
%TRIANGLESTANDARDFORM_REF puts a triangle in standard form with respect to a
%given edge for refienment purposes
%tri: 3 indices defining a triangle
%edge: 2 indices defining an edge
%If the edge is not contained in the triangle, the triangle in unchanged.


tri_out = tri;
fledge = flip(edge);


for i = 1:3
    
    t = circshift(tri_out, i);
    
    if sum(t(2:3) == edge) == 2
        
        tri_out = circshift(t, 1); %Makes sure that the edge of interest is between vertices 3 and 1
        
    elseif sum(t(2:3) == fledge) == 2
        
        tri_out = circshift(t, 1); %Makes sure that the edge of interest is between vertices 3 and 1
        
    end
    
end





end

