%
% This code accompanies the paper:
%
% "Partial Functional Correspondence"
% Rodola, Cosmo, Bronstein, Torsello, Cremers
% Computer Graphics Forum 2016
%
% Please cite the paper above if you use this code in your research.
%
% Written by Emanuele Rodola and Luca Cosmo
%
function [C, est_rank, final_cost] = optimize_F(M, N, G, F, v, C_init, mu1, mu2)

fprintf('Optimizing C...\n');

k = size(M.BASIS.basis, 2);

% COMPUTE W MASK
est_rank = sum(N.BASIS.eigenvalues-max(M.BASIS.eigenvalues)<0);
W = zeros(k);
for i=1:k
    for j=1:k
        slope = est_rank/k;
        direction = [1 slope];
        direction = direction./norm(direction);
        W(i,j) = exp(-0.03*sqrt(i.^2 + j.^2))*norm(cross([direction 0], [i,j, 0]-[1 1 0]));        
    end
end
% figure, imagesc(W);
d=zeros(1,k);
d(1:est_rank)=1;

D = repmat(d,k,1);

fprintf('Estimated rank: %d\n', est_rank);

if isempty(C_init)
    C_init = (max(max(W))-W)./max(max(W));
end

% 

A = N.BASIS.basis'*N.SHAPE.A*F;
B = M.BASIS.basis'*M.SHAPE.A*diag(sparse(v))*G;

manifold = euclideanfactory(k,k);
problem = {};

problem.M = manifold;

problem.cost = @(C) (...
    sum(sum((C*A-B).^2).^0.5) + ...
    mu1 * norm(C.*W,'fro')^2 + ...
    mu2 * (norm(C'*C,'fro')^2 - sum(diag(C'*C).^2) + sum((diag(C'*C) - d').^2) ));

problem.egrad = @(C) (...
    norm_21_gradient(C,A,B) + ...
    mu1 * 2 * C.*W.*W + ...
    mu2 * 4*(C*C'*C - C.*D ));

% figure, checkgradient(problem)

options.maxiter = 1e4;
options.tolgradnorm = 1e-06;
options.minstepsize = 1e-06;
options.verbosity = 1;

[C, final_cost, info, ~] = conjugategradient(problem, C_init, options);
fprintf('C-step, iterations: %d, cost: %f\n', info(end).iter, final_cost);

end
