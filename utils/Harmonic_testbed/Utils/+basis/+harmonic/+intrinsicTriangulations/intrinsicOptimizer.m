function [Basis, radii1] = intrinsicOptimizer(shape_number, S1, Basis2, landmarks1, landmarks2, NumControlLandmarks, radii1)
%CENTRALFEOPTIMIZER produces the optimized basis 

radii1 = radii1(:);

funObj = @(v) intrinsicEnergy(shape_number, S1, Basis2, landmarks1, landmarks2, NumControlLandmarks, v);
funProj = @(F) F;


options.verbose = 1;
options.maxIter = 500;
% options.optTol = 1e-15; % 1e-5 default.
% options.progTol = 1e-15; % 1e-9 default.
options.numDiff = 1; %Finite Differences

% size(radii1)
radii1 = minConf_PQN( funObj, radii1, funProj, options);
% size(radii1)

Basis = computeIntrinsicBasis( shape_number, S1.nv, landmarks1(1:end-NumControlLandmarks), radii1 );









end