function surface_out = insertVertexInTriangle(surface, triangle_ind)
%INSERTVERTEXINTRIANGLE inserts a vertex in the center of a triangle,
%updating the list of vertices and triangles
%surface: structure with fields TRIV, X, Y, Z and VERT
%triangle_ind: index of the triangle to be refined

surface_out = surface;


vertices_tri = surface.TRIV(triangle_ind, :);

center = ( surface.VERT(vertices_tri(1),:) + surface.VERT(vertices_tri(2),:) + surface.VERT(vertices_tri(3),:) )/3;


surface_out.VERT = [surface_out.VERT; center]; % adds new vertex
center_ind = size(surface_out.VERT, 1);


surface_out.TRIV(triangle_ind,:) = [vertices_tri(1) vertices_tri(2) center_ind]; %Replaces old triangle with first new triangle. This preserves (part of) the triangle indexing.

surface_out.TRIV = [surface_out.TRIV;... %Places the other two new triangles at the end of the list.
                    vertices_tri(2) vertices_tri(3) center_ind;...
                    vertices_tri(3) vertices_tri(1) center_ind];




surface_out.X = surface_out.VERT(:,1); % Updates the separate coordinates.
surface_out.Y = surface_out.VERT(:,2);
surface_out.Z = surface_out.VERT(:,3);


end

