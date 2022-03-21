function funcs_out = BoundaryReinsertBatch(funcs, boundary_idx, boundary_values)
%Reinserts boundary values into arrays of functions.

funcs_out = zeros(size(funcs,1) + length(boundary_idx), size(funcs,2));

for i = 1:size(funcs,2)
    
    funcs_out(:,i) = BoundaryReinsert(funcs(:,i), boundary_idx, boundary_values);   
    

end





end