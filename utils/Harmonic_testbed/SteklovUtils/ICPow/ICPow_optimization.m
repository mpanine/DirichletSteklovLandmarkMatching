function [opt_p2p, opt_Basis2] = ICPow_optimization(Basis1, Basis2, numiter, numskip)
%ICPow_optimization optimizes Basis2 via an ICP-like power method.
%Outputs the optimized version of Basis2 and the optimized p2p map S2 -> S1.

%Basis1: harmonic basis on first shape.
%Basis2: harmonic basis on second shape. Powers are applied here.
%numiter: number of iterations of the algorithm. (Optional argument)
%numskip: used to speed up the execution by skipping points. (Optional argument)

if nargin<4 %Default values
    numskip = 1;
end
 
if nargin<3 %Default values
    numiter = 5;
end


opt_p2p = annquery(Basis1', Basis2', 1); 
opt_Basis2 = Basis2;


for i = 1:numiter
    
    [opt_p2p, opt_Basis2] = ICPow_iteration(Basis1, opt_Basis2, opt_p2p, numskip);   
    
end




end