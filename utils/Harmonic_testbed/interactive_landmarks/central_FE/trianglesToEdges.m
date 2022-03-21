function edges = trianglesToEdges(triangles)
%TRIANGLESTOEDGES extracts all the distinct non-oriented edges from a list
%of triangles.
%triangles: list of triangles (size = (ntriangles, 3) )

edges = zeros(3*size(triangles,1), 2);

ed = 1;

for t = 1:size(triangles,1)
    
    edges(ed,:) = [triangles(t,1) triangles(t,2)];
    ed = ed + 1;
    edges(ed,:) = [triangles(t,2) triangles(t,3)];
    ed = ed + 1;
    edges(ed,:) = [triangles(t,3) triangles(t,1)];
    ed = ed + 1;

end

edges = sort(edges,2);
edges = unique(edges,'rows');


end