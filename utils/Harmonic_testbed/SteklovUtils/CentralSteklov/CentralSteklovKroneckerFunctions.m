function [H, ND] = CentralSteklovKroneckerFunctions(Shape_refined)
%CentralSteklovKroneckerFunctions: produces the expansions of the harmonic
%functions with Kroncker delta boundary conditions in the Steklov basis.
%Also produces the expansions for the relevant 

num_bounds = length(Shape_refined.STEKLOV.evecs_bound);

H = zeros(size(Shape_refined.STEKLOV.evals, 1), num_bounds);
ND = zeros(size(Shape_refined.STEKLOV.evals, 1), num_bounds);

for i = 1:num_bounds %Cycle over all Kronecker delta BC
    
    BC = zeros(1,num_bounds);
    BC(i) = 1;    
    
    for j = 1:num_bounds %Cycle over all boundaries
        
%         size(ones(1, size(Shape_refined.STEKLOV.evecs_bound{j}, 1) ) * Shape_refined.STEKLOV.MassBoundary{j} * Shape_refined.STEKLOV.evecs_bound{j})
        
%        H(:,i) = H(:,i) + BC(j) * ones(size(Shape_refined.STEKLOV.evecs_bound{j}))' * Shape_refined.STEKLOV.MassBoundary{j} * Shape_refined.STEKLOV.evecs_bound{j} ;       
        H(:,i) = H(:,i) +  BC(j) * transpose(ones(1, size(Shape_refined.STEKLOV.evecs_bound{j}, 1) ) *...
                                    Shape_refined.STEKLOV.MassBoundary{j} * Shape_refined.STEKLOV.evecs_bound{j} );  
                                
        ND(:,i) = ND(:,i) +  BC(j) * transpose(ones(1, size(Shape_refined.STEKLOV.evecs_bound{j}, 1) ) *...
                                    Shape_refined.STEKLOV.MassBoundary{j} * Shape_refined.STEKLOV.evecs_bound{j} * Shape_refined.STEKLOV.evals); %Normal derivatives via the eigenvalues
                                
    end  
    
    
    

    
    
end


end