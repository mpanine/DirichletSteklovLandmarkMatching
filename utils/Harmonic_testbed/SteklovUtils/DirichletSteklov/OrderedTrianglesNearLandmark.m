function tri_list_new = OrderedTrianglesNearLandmark(surface, landmark )
%Finds the triangles adjacent to a landmark and orders them in a consitent
%manner. This is used to automatically order the new vertices in the
%landmark refinement.

%surface: SHAPE.surface
%index of the landmark under consideration

tri_list = findTrianglesContainingVertex(surface, landmark );


tri_list_placed = zeros(size( tri_list )); %Tracks whether the triangle has been placed.
tri_list_placed(1) = 1;


tri_list_new = tri_list(1); % The placement of the first triangle does not matter.



while sum(tri_list_placed) ~= length(tri_list)
    
    for i = 1:length(tri_list)
        if ~tri_list_placed(i) % The triangle is not already placed.
            
            %Check if triangle i can be placed to the left of the list.
            if checkTriangleAdjacence(surface.TRIV, tri_list_new(1), tri_list(i) ) %Is it adjacent to the left end of the list?          
                
                    if checkLeftRightTriangleAdjacence(surface.TRIV, tri_list_new(1), tri_list(i), landmark) == -1
                        tri_list_new = [tri_list(i); tri_list_new];
                        tri_list_placed(i) = 1;
                    end

            end

            if ~tri_list_placed(i) % Must be checked again to avoid double placing the last triangle
                %Check if triangle i can be placed to the right of the list.
                if checkTriangleAdjacence(surface.TRIV, tri_list_new(end), tri_list(i) ) %Is it adjacent to the right end of the list?          

                        if checkLeftRightTriangleAdjacence(surface.TRIV, tri_list_new(end), tri_list(i), landmark) == 1
                            tri_list_new = [tri_list_new; tri_list(i)];
                            tri_list_placed(i) = 1;
                        end

                end
            end
        
        end
    end
        
        
    
    
    
    
    
    
    
end


end