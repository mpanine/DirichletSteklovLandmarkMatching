function vertices = trianglesToVertices(triangles)
%TRIANGLESTOVERTICES extracts all the distinct vertices from a list
%of triangles.
%triangles: list of triangles (size = (ntriangles, 3) )

vertices = unique( triangles );

end