function [ metrics ] = map_dirichlet_energy(Src,Tar,pF)
%COMPUTES THE DIRICHLET ENERGY (the energy minimized by harmonic maps) of a
%point-to-point map between two shapes
%Src: source shape
%Tar: target shape
%pF: point-to-point map Tar -> Src (a vector of length numvertex(Tar) )

%NOTE: the statistical quantitites and the associated terminology do not
%make sense for this energy. Namely:
%mean: really the value of the Dirichlet energy
%standardDeviation: a sort of equivalent of standard deviation, assuming
    %that mean is an expectation value in the quantum-mechanical sense (without normalization)
%worstCase: does not make much sense, either.



    assert(isfield(Tar, 'Gamma'),...
    'The Source shape needs to have a Gamma field to compute the distortion metrics. Set computeGeoDist to true when creating the shape.');
    assert(isfield(Tar, 'Gamma'),...
    'The Target shape needs to have a Gamma field to compute the distortion metrics. Set computeGeoDist to true when creating the shape.');


    D = Tar.W .* Src.Gamma(pF,pF).^2 ;
    
    
    metrics.mean = full( sum(sum( D )) );
%     metrics.values = D;
    metrics.values = [];
%     % TODO understand how to compute std on this
%     metrics.standardDeviation =  sqrt(abs(   full(sum(sum(  Tar.W^2 .* Src.Gamma(pF,pF).^2  )))  -   metrics.mean^2   ));    
%     metrics.worstCase = full( max(max( D )) );



end

