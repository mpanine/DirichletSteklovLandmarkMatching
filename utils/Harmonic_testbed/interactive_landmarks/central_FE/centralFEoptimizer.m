function radii1 = centralFEoptimizer(S1, Basis2, landmarks1, landmarks2, NumControlLandmarks, radii1)
%CENTRALFEOPTIMIZER optimizes the radii of the central elements of S1 w.r.t
%the values of Basis2 at the control landmarks





funObj = @(v) centralFEenergy(S1, Basis2, landmarks1, landmarks2, NumControlLandmarks, v);        
funProj = @(F) F;


options.verbose = 1;
options.maxIter = 50;
% options.optTol = 1e-15; % 1e-5 default.
% options.progTol = 1e-15; % 1e-9 default.
options.numDiff = 1; %Finite Differences


radii1 = minConf_PQN( funObj, radii1, funProj, options);










end