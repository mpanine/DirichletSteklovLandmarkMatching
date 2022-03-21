function out = combineProductFmapCellToMatrix(Fmap_cell_1, Fmap_cell_2, iftranspose )
%Combines the products of two cell arrays of matrices into one matrix. Useful for illustrating
%Dirichlet-Steklov Fmaps on multiple boundaries.


num_bounds = length(Fmap_cell_1);

out = [];

for i = 1:num_bounds
    
    if iftranspose    
        out = [out; Fmap_cell_1{i}' * Fmap_cell_2{i}];      
    else
        out = [out; Fmap_cell_1{i} * Fmap_cell_2{i}];          
    end
    
    
end








end