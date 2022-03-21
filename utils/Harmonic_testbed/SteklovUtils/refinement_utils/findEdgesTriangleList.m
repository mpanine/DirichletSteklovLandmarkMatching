function edge_list = findEdgesTriangleList(surface, tri_list)
%FINDEDGESTRIANGLELIST finds the edges of a list of triangles. Removes
%duplicates
%surface: structure containing field TRIV
%tri_list: contains list of triangles

edge_list = [];


for t = 1:length(tri_list)
    
    tri = surface.TRIV(tri_list(t), :);
    edge_list = [edge_list; sort(tri(1:2)); sort(tri(2:3)); sort([tri(3) tri(1)])];
      
    
end

edge_list = unique( edge_list, 'rows' );



end

