function [next_p2p, next_Basis2] = ICPow_iteration(Basis1, Basis2, current_p2p, numskip)
%ICPow_iteration produces a single iteration of the ICPow iteration.
%The goal is to obtain powers that align Basis1^powers to Basis2
%This code is executed iteratively in an ICP-like fashion. We term this
%ICPow.

%Basis1: harmonic basis on first shape.
%Basis2: harmonic basis on second shape. Powers are applied here.
%current_p2p: point-to-point map at the beginning of the current iteration. Goes Shape2 -> Shape1

if nargin<4
    numskip = 1; %Default no point skip.
end

if nargin<3   
   current_p2p = annquery(Basis1', Basis2', 1);    
end


powers = ones(1, size(Basis1, 2));   

funObj = @(v) ICPowDistance(Basis1, Basis2, v, current_p2p, numskip);
funProj = @(F) F;

options.verbose = 0;
options.maxIter = 500;
options.optTol = 1e-15; % 1e-5 default.
options.progTol = 1e-15; % 1e-9 default.


% [InitialDistance, InitGrad] = controlLandmarkDistance(Basis1, Basis2, landmarks1, landmarks2, NumControlLandmarks, powers);
% InitialDistance
% InitGrad;

next_powers = minConf_PQN( funObj, powers', funProj, options);
next_powers = next_powers';



next_Basis2 = Basis2.^next_powers;
next_p2p = annquery(Basis1', next_Basis2', 1);


end