function reorder = boundaryReorder(boundary_edges)
% reorders the boundary in order to make it possible to plot the functions
% on it more or less smoothly. This should split the boundary into circles
% and segments (assuming that everything works as intended).

%reorder: cell array with boundary components as different cells

numedges = size(boundary_edges, 1);
ifplaced = zeros(numedges,1);


good_edge_found = false; % Used in finding a starting edge

for tt = 1:numedges
    
    if findNumEdgeNeighborsOfEdge(boundary_edges, tt) < 3        
        reorder{1} = boundary_edges(tt,:); % The first edge starts the order.
        ifplaced(tt) = 1;
        good_edge_found = true;
        break;
    end
        
end

if ~good_edge_found    
    warning('No good starting edge found. Starting wherever.')
    reorder{1} = boundary_edges(1,:); % The first edge starts the order
    ifplaced(1) = 1;     
end

current_component = 1; %Index of the current connected component.
num_placed_edges = 1;


while true % The cycle runs until all edges are placed. (see break below)
    
    current_placed_edges = num_placed_edges;
    
    for i = 1:numedges % Cycle over all edges
        if ifplaced(i) ~= 1 % Unplaced edges only
            
            if boundary_edges(i,1) == reorder{current_component}(end)
                
                reorder{current_component} = [reorder{current_component} boundary_edges(i,2)];
                ifplaced(i) = 1;
                num_placed_edges = num_placed_edges + 1;
                
            elseif boundary_edges(i,2) == reorder{current_component}(end)
                
                reorder{current_component} = [reorder{current_component} boundary_edges(i,1)];
                ifplaced(i) = 1;
                num_placed_edges = num_placed_edges + 1;
                
            elseif boundary_edges(i,1) == reorder{current_component}(1)
                
                reorder{current_component} = [boundary_edges(i,2) reorder{current_component}];
                ifplaced(i) = 1;
                num_placed_edges = num_placed_edges + 1;
                
            elseif boundary_edges(i,2) == reorder{current_component}(1)
                
                reorder{current_component} = [boundary_edges(i,1) reorder{current_component}];
                ifplaced(i) = 1;
                num_placed_edges = num_placed_edges + 1;
                
            end
            
        end
    end
    
    if num_placed_edges == current_placed_edges %A full cycle over unplaced edges has failed to place new edges.
        
        if num_placed_edges == numedges % We are done
            
            break;
            
        else % We need a new component
            
            current_component = current_component + 1;
            thingy = find(ifplaced == 0);
            
            good_edge_found = false;
            
            for tt = 1:length(thingy)
                
                    if findNumEdgeNeighborsOfEdge(boundary_edges, thingy(tt)) < 3        
                        reorder{current_component} = boundary_edges(thingy(tt),:); % The first edge starts the order.
                        ifplaced(thingy(tt)) = 1;
                        num_placed_edges = num_placed_edges + 1;
                        good_edge_found = true;
                        break;
                    end
                    
           
            end
            
            if ~good_edge_found
                warning('No good starting edge found. Starting wherever.')
                reorder{current_component} = boundary_edges(thingy(1),:); % The first edge starts the order
                ifplaced(thingy(1)) = 1;
                num_placed_edges = num_placed_edges + 1;
            end
            
            
            
        end
        
    end

end

    


end


function numed = findNumEdgeNeighborsOfEdge(boundary_edges, edge_ind)
%Finds the number of edges neighboring a given edge.

other_edges = boundary_edges;
other_edges(edge_ind) = [];

numed = sum( ismember(other_edges, edge_ind) );

end

