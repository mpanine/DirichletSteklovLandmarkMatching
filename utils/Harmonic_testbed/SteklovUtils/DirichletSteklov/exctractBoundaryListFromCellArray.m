function B = exctractBoundaryListFromCellArray(boundaries, component_list)
%boundaries: cell array of lists of vertices (boundary components)
%component_list: list of components to combine

B = [];

for i = 1:length(component_list)
    
   B = [B; boundaries{i}]; 
    
    
end








end