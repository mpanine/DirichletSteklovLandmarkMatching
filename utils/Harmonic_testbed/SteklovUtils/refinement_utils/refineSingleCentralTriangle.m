function [shape, normalTriangles] = refineSingleCentralTriangle(shape, triangle, center, radius, numsubdivisions)
%refineSingleCentralTriangle refines a triangle containing the vertex
%center in order to create a small "circle" of radius "radius" around the
%center.
%shape: struncture with field surface
%triangle: index of the triangle (refers to shape.surface.TRIV)
%center: vertex used as the center of the circle being built
%radius: radius of the circle being built
%numsubdivisions: number of subdivisions of the circle being built. The
%more, the smoother. The minimum is 1 (no subdivision).

if numsubdivisions < 1 %There should be at least one part
    numsubdivisions = 1;
end

nv = size(shape.surface.VERT,1);


tri = shape.surface.TRIV(triangle,:);

for s = 1:3 % This loop converts the triangle to standard form   
    if tri(1) == center
       break 
    else
       tri = circshift(tri,1);
    end
    
    if s == 3
        error('The center vertex is not contained in the triangle')
    end    
end

% Normalized vectors from the center to the existing vertices
V1 = shape.surface.VERT(tri(2),:) - shape.surface.VERT(tri(1),:);
V1 = V1/norm(V1);
V2 = shape.surface.VERT(tri(3),:) - shape.surface.VERT(tri(1),:);
V2 = V2/norm(V2);

V3 = V2 - V1*(V1*V2'); % Normalized vector orthogonal to V1 and coplanar to V1 and V2.
V3 = V3/norm(V3);


%% Building the new vertices

maxangle = acos(V1*V2');
angles = linspace(0, maxangle, numsubdivisions + 1);%Angles going from V1 to V2

circleVertices = zeros(numsubdivisions + 1, 3);

for i = 1:numsubdivisions + 1    
   circleVertices(i,:) = shape.surface.VERT(center,:) + radius*( cos(angles(i))*V1 + sin(angles(i))*V3 );    
end


%Now to build the vertices on the opposite edge of the original triangle.

l1 = norm( shape.surface.VERT(tri(2),:) - shape.surface.VERT(tri(1),:) ); %Quantities for Rmax
a = (shape.surface.VERT(tri(3),:) - shape.surface.VERT(tri(1),:)) * V1'; 
b = (shape.surface.VERT(tri(3),:) - shape.surface.VERT(tri(1),:)) * V3'; 

Rmax = @(t) l1/( cos(t) - (a - l1)/b * sin(t) );

outerVertices = zeros(numsubdivisions - 1,3);

for i = 1:size(outerVertices,1)
   
    outerVertices(i,:) = shape.surface.VERT(center,:) + Rmax(angles(i+1))*( cos(angles(i+1))*V1 + sin(angles(i+1))*V3 );
    
end


%% Find the other triangle connected to edge opposite to center (ASSUMES MANIFOLD: AT MOST ONE SUCH TRIANGLE EXISTS)

[tri_indices, ~] = findTrianglesContainingEdge(shape.surface, [tri(2) tri(3)]);

second_triangle = [];


for i = 1:length(tri_indices)    
   if tri_indices(i) ~= triangle
       second_triangle = tri_indices(i);    
   end    
end

if ~isempty(second_triangle)
   
    second_tri = shape.surface.TRIV(second_triangle,:);
    
    for i=1:3 % Find the vertex not part of the common edge, call it co_center
       if (second_tri(i) ~= tri(2) ) && (second_tri(i) ~= tri(3))
           co_center = second_tri(i);
       end
    end
    
    for s = 1:3 % This loop converts the triangle to standard form   
        if second_tri(1) == co_center
           break 
        else
           second_tri = circshift(second_tri,1);
        end
    
        if s == 3
            error('The co_center vertex is not contained in the second triangle. This is super bad.')
        end    
    end
    
end






 



%% Indexing the new vertices
% Circle   -- Must make sure to not create multiple copies of the same vertex

threshold = radius/1000; %Threshold for vertex coincidence.

% Check if the circle vertices are duplicates of already existing vertices

[D1, P1] = min( sqrt( sum(( shape.surface.VERT - circleVertices(1,:)).^2  , 2) ) ); % Finding the closest vertex.

if D1 < threshold
    ifFirstDuplicate = true;
    firstIndex = P1;
    circleVertices(1,:) = [];    
else    
    ifFirstDuplicate = false;
end


[D2, P2] = min( sqrt( sum(( shape.surface.VERT - circleVertices(end,:)).^2  , 2) ) ); % Finding the closest vertex.

if D2 < threshold
    ifLastDuplicate = true;
    lastIndex = P2;
    circleVertices(end,:) = [];    
else
    ifLastDuplicate = false;
end

circleVertexInds = [nv + 1:nv + size(circleVertices,1)]; % This only includes the new vertices (duplicates have been removed)

if ifFirstDuplicate
    circleVertexInds = [P1 circleVertexInds]; 
end

if ifLastDuplicate
    circleVertexInds = [circleVertexInds P2]; 
end




% Outer indices
outerVertexInds = [nv + size(circleVertices,1) + 1: nv + size(circleVertices,1) + size(outerVertices,1)];
outerVertexInds = [tri(2) outerVertexInds tri(3)]; %Include the indices of the original vertices





%% Building the new triangles in triangle adjacent to center



% Inner circle

innerTriangles = zeros(numsubdivisions, 3);

for i = 1:numsubdivisions    
    innerTriangles(i,:) = [center circleVertexInds(i) circleVertexInds(i+1)];   
end


%Outer Circle

outerTriangles = zeros(numsubdivisions*2, 3);
normalTriangles = zeros(numsubdivisions, 3); % Triangles used in the computation of the normal derivatives.

for i = 1:numsubdivisions
    
   outerTriangles(2*i-1,:) = [circleVertexInds(i) outerVertexInds(i) outerVertexInds(i+1)]; % Two outer vertices.
%    outerTriangles(2*i,:)   = [circleVertexInds(i) outerVertexInds(i+1) circleVertexInds(i+1)];  % Two circle vertices.
   outerTriangles(2*i,:)   = [circleVertexInds(i+1) circleVertexInds(i) outerVertexInds(i+1)];  % Two circle vertices: can be used for normal derivative computation. 
   normalTriangles(i,:) = [circleVertexInds(i+1) circleVertexInds(i) outerVertexInds(i+1)];
   
end


%% Build new triangles in triangle adjacent to the edge opposite to the center

if ~isempty(second_tri)
    
    exteriorTriangles = zeros(numsubdivisions, 3);
    
    for i = 1:numsubdivisions    
        exteriorTriangles(i,:) = [outerVertexInds(i) co_center outerVertexInds(i+1)];   
    end
    
end

%% Update the shape structure

%Expand the list of vertices

shape.surface.VERT = [shape.surface.VERT; circleVertices; outerVertices];
shape.surface.X = shape.surface.VERT(:,1);
shape.surface.Y = shape.surface.VERT(:,2);
shape.surface.Z = shape.surface.VERT(:,3);

%Update the triangle list

shape.surface.TRIV(triangle,:) = innerTriangles(1,:);
shape.surface.TRIV = [shape.surface.TRIV; innerTriangles(2:end,:); outerTriangles];


%Insert the exterior triangles, if any

if ~isempty(second_tri)
    
    shape.surface.TRIV(second_triangle,:) = exteriorTriangles(1,:);
    shape.surface.TRIV = [shape.surface.TRIV; exteriorTriangles(2:end,:)];
    
end
























end