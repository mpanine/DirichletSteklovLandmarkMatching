function func_rebuilt = BoundaryReinsert(func_cut, boundary_idx, boundary_values)
%Reinserts the boundary elements into a function without boundary elements
%func_cut: function of non-boundary vertices of a mesh
%boundary_idx: indices of the boundary vertices
%boundary_values: values to be inserted at the boundary (boundary conditions, usually)

func_rebuilt = zeros(1,length(func_cut) + length(boundary_idx));

current_ind = 0;

for i = 1:length(func_rebuilt)
    
   if sum(i == boundary_idx) == 0  % i is not a boundary point
       
       current_ind = current_ind + 1;
       func_rebuilt(i) = func_cut(current_ind);
        
   else % i is a boundary point
        
%        func_rebuilt(i) = boundary_values( find(boundary_idx == i) );
       func_rebuilt(i) = boundary_values( boundary_idx == i );
    
   end
    
    
end

end
