function a = checkTriangleAdjacence(TRIV, t1, t2)
%Checks whether triangles t1 and t2 of TRIV are adjacent.

%TRIV:list of triangles
%t1: index of first triangle
%t2: index of second triangle


T1 = TRIV(t1,:);
T2 = TRIV(t2,:);


num_common_vert = sum( ismember(T1,T2) );

if num_common_vert == 2
    a = true;
else
    a = false; % Identical triangles are not adjacent for our purposes!
end







end