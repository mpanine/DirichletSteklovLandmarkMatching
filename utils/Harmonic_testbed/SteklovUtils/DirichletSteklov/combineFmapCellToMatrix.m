function out = combineFmapCellToMatrix(Fmap_cell)
%Combines a cell array of matrices into one matrix. Useful for illustrating
%Dirichlet-Steklov Fmaps on multiple boundaries.


num_bounds = length(Fmap_cell);

out = [];

for i = 1:num_bounds
    
    out = [out; Fmap_cell{i}];  
    
    
end








end