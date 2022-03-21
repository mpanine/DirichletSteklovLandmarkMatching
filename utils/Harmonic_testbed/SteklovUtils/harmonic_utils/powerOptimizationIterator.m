function Current_Basis_2 = powerOptimizationIterator(Basis1, Basis2, numlandmarks, sigma, num_iterations, options)
%POWEROPTIMIZATIONITERATOR obtains the optimal embedding after power
%optimization iteration




funObj = @(v) kernelDistance(Basis1, Basis2, sigma, v);
        
funProj = @(F) F;
        

powers = ones( numlandmarks , 1);
powers = minConf_PQN( funObj, powers, funProj, options);
powers = powers';


Current_Basis_2 = Basis2.^powers;

for iter = 1:num_iterations
    
    powers_iter = 0.01*ones( numlandmarks , 1); %
    funObj_iter = @(v) iterativeKernelDistance(Basis1, Current_Basis_2, sigma, v);
    powers_iter = minConf_PQN( funObj_iter, powers_iter, funProj, options);
    
    Current_Basis_2 = 0.5 * Current_Basis_2 + 0.5 * Current_Basis_2.^(powers_iter' + 1) ;
    
end






end

