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
function v = optimize_v(M,N,G,F,C,mu1,mu2,opt)
%M = Src
%N = Tar
%G = descr_Src
%F = descr_Tar
% M.A = M.Ae
    fprintf('Optimizing v...\n');
    
    k = size(C,1);

    A = C*N.BASIS.basis'*N.SHAPE.A*F;
    B = M.BASIS.basis'*M.SHAPE.A;

    areas = full(diag(M.SHAPE.A));
    target_area = full( sum(diag(N.SHAPE.A)) );

    x0 = M.BASIS.basis*(C*(N.BASIS.basis'*(N.SHAPE.A*ones(N.SHAPE.nv,1))));
  
    manifold = euclideanfactory(M.SHAPE.nv,1);

    problem = {};
    problem.M = manifold;

    vfunc = mumford_shah(M.SHAPE.surface.VERT, M.SHAPE.surface.TRIV, M.SHAPE.A);
    problem.cost =  @(v) vfunc.cost(A, B, G, v, ones(size(v)), target_area, areas, mu1, mu2, opt);
    problem.egrad = @(v) vfunc.grad(A, B, G, v, ones(size(v)), target_area, areas, mu1, mu2, opt);

%     figure, checkgradient(problem);

    options.maxiter = 2e2;%5e2;
    options.tolgradnorm = 1e-6;
    options.minstepsize = 1e-6;
    options.verbosity = 1;

    [v, cost, info, ~] = conjugategradient(problem, x0, options);
	fprintf('v-step, iterations: %d, cost: %f\n', info(end).iter, cost);
end
