function [boundary_new, boundary_edges_old] = findCentralLandmarkBoundary(Shape, segment, boundary_old)
% Finds the boundary of a segment of a shape for central landmark refinement.
%Shape: structure containing fields .surface.TRIV and .surface.VERT
%segment: list of vertices of Shape that form the desired segment.

%boundary_old: boundary indices in the origial indices 
%boundary_new: boundary indices in the indices relative to the order of "segment"
%boundary_edges_old: defines the boundary connectivity



numtri = size( Shape.SHAPE.surface.TRIV, 1);

all_edges = findEdgesTriangleList(Shape.SHAPE.surface, 1:numtri);

pre_boundary_edges = sum( ismember( all_edges, boundary_old ), 2); 
boundary_edges_old = all_edges( pre_boundary_edges == 2 , : ); % Both vertices of edge in boundary.

% size(boundary_edges_old)

% pre_segment_edges = sum( ismember( all_edges, segment ), 2); 
% boundary_cross_edges = all_edges( pre_segment_edges == 1, :); % Only one vertex of the edge is in the segment.
% 
% 
% 
% boundary_cross_vertices = boundary_cross_edges(:);
% boundary_cross_vertices = unique(boundary_cross_vertices); % All vertices on either side of a boundary.
% 
% 
% boundary_old = boundary_cross_vertices( ismember(boundary_cross_vertices, segment) );

aa = 1:length(segment);
boundary_new = aa( ismember( segment, boundary_old) ); 
    







%% The boundary must not contain triangles

% Boundary edges must be part of at most one triangle in segment.
% To speed up the search: find boundary vertices part of more than two boundary edges



preBadTri = sum( ismember(Shape.SHAPE.surface.TRIV, boundary_old) , 2);
BadTri = find(preBadTri == 3); % List of triangles fully contained in the boundary.

if ~isempty(BadTri)
    
    thingy = [1 2 3 1];
    problematic_edges = [];%indices of the problematic edges in the boundary.

    for i = 1:length(BadTri)
    
        for j = 1:3
    
            edge = [Shape.SHAPE.surface.TRIV(BadTri(i), j) Shape.SHAPE.surface.TRIV(BadTri(i), thingy(j + 1) )];
            [tri_indices, ~] = findTrianglesContainingEdge(Shape.SHAPE.surface, edge);
            
            num_tri_in_segment = 0;
            
            for k = 1:length(tri_indices)                
                tri = Shape.SHAPE.surface.TRIV(tri_indices(k),:);
                num_vert_in_segment = sum( ismember(tri,segment) );
                
                if num_vert_in_segment == 3
                    num_tri_in_segment = num_tri_in_segment + 1;
                end                
            end
            
            if num_tri_in_segment > 1 % This edge is between two triangles in segment.
                
                bad_edge = find(ismember(boundary_edges_old, edge,'rows'));
                if isempty(bad_edge)
                    bad_edge = find(ismember(boundary_edges_old, circshift(edge,1),'rows'));
                end
                                
                problematic_edges = [problematic_edges bad_edge];
            end
            
    
        end
    end
        
    boundary_edges_old(problematic_edges,:) = []; % Remove the problematic edges.
    
end




%BACKUP VERISON 2  -- Decent

% potentially_bad_edges = [];
% 
% for i = 1:length(boundary_old)
%     
%     numedges = sum(sum( boundary_edges_old == boundary_old(i) ));
%     
%     if numedges > 2       
%         [edge_list,~] = find(boundary_edges_old == boundary_old(i));
%         potentially_bad_edges = [potentially_bad_edges; edge_list];
%     end
%     
% end
% 
% problematic_edges = [];
% 
% if ~isempty(potentially_bad_edges)
% 
%     for i = 1:length(potentially_bad_edges) 
% 
%         [tri_indices, ~] = findTrianglesContainingEdge(Shape.surface, boundary_edges_old(potentially_bad_edges(i),:));
% 
%         num_tri_in_segment = 0;
% 
%         for j = 1:length(tri_indices)
% 
%             tri = Shape.surface.TRIV(tri_indices(j),:);
%             num_vert_in_segment = sum( ismember(tri,segment) );
% 
%             if num_vert_in_segment == 3
%                 num_tri_in_segment = num_tri_in_segment + 1;
%             end
% 
%         end
% 
%         if num_tri_in_segment > 1
%             problematic_edges = [problematic_edges potentially_bad_edges(i)];
%         end
% 
% 
%     end
% 
% end
% 
% boundary_edges_old(problematic_edges,:) = []; % Remove the problematic edges.

%BACKUP SLOW VERSION

% problematic_edges = [];
% for i = 1:size(boundary_edges_old,1)  %NOTE: THIS IS SLOW. FIND A WAY TO SPEED THIS UP.
% 
%     [tri_indices, ~] = findTrianglesContainingEdge(Shape.surface, boundary_edges_old(i,:));
%     
%     num_tri_in_segment = 0;
%     
%     for j = 1:length(tri_indices)
%         
%         tri = Shape.surface.TRIV(tri_indices(j),:);
%         num_vert_in_segment = sum( ismember(tri,segment) );
%         
%         if num_vert_in_segment == 3
%             num_tri_in_segment = num_tri_in_segment + 1;
%         end
%         
%     end
%     
%     if num_tri_in_segment > 1
%         problematic_edges = [problematic_edges i];
%     end
%     
% 
% end
% 
% boundary_edges_old(problematic_edges,:) = []; % Remove the problematic edges.



end