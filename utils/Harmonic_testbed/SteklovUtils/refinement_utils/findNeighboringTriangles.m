function findNeighboringTriangles(surface, triangle_ind)
%FINDNEIGHBRINGTIRANGLES outputs the indices of the triangles adjacent to a
%target triangle. Adjacence is only considered across edges.
%surface: structure with fields TRIV, X, Y, Z and VERT
%triangle_ind: index of the triangle whose neighbors we are trying to find.


center_tri = surface.TRIV(triangle_ind, :); %Vertices of the central triangle

center_edges = [center_tri(1) center_tri(2); center_tri(2) center_tri(3); center_tri(3) center_tri(1)]; %Edges of center triangle
    













end