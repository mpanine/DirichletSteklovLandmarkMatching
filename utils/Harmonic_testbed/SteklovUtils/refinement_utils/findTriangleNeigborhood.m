function [tri_list, all_vertices] = findTriangleNeigborhood(surface, center, num_circles)
%FINDTRIANGLENEIGHBORHOOD builds a list of triangles that are num_circles
%neighbors of a vertex center
%shape: structure with field surface
%center: index of the center vertex
%num_circles: number of "circles of triangles" to be found. Always does at
%least one circle

tri_list = findTrianglesContainingVertex(surface, center); %first circle


for n = 2:num_circles %Other circles

    current_vertices = surface.TRIV(tri_list, :);
    current_vertices = current_vertices(:);
    current_vertices = unique(current_vertices);
    
    for v = 1:length(current_vertices)
        
        tri_list = [tri_list ; findTrianglesContainingVertex(surface, current_vertices(v)) ];
    
    end


end

tri_list = unique(tri_list);

all_vertices = surface.TRIV(tri_list, :);
all_vertices = all_vertices(:);
all_vertices = unique(all_vertices);



end