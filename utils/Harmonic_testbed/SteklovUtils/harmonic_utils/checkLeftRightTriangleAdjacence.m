function direction = checkLeftRightTriangleAdjacence(TRIV,t1,t2,center)
%Outputs the side of t1 to which t2 is adjacent. THE TRIANGLES MUST
%ALREADY BE ADJACENT!
%Both triangles MUST contain the vertex "center" and the direction is determined
%relative to it and the common orientation of the triangles.

% left = -1
% right = 1


%% Put the triangles in standard form

T1 = TRIV(t1,:);

for s = 1:3 % This loop converts the triangle to standard form   
    if T1(1) == center
       break 
    else
       T1 = circshift(T1,1);
    end
    
    if s == 3
        error('The center vertex is not contained in the triangle')
    end    
end


T2 = TRIV(t2,:);

for s = 1:3 % This loop converts the triangle to standard form   
    if T2(1) == center
       break 
    else
       T2 = circshift(T2,1);
    end
    
    if s == 3
        error('The center vertex is not contained in the triangle')
    end    
end


if T1(3) == T2(2)
    direction = 1;
elseif T1(2) == T2(3)
    direction = -1;
else
    error('The triangles must share an edge containing the center!')
end









end