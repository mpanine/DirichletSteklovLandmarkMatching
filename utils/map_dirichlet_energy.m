function [ metrics ] = map_dirichlet_energy(Src,Tar,Tar_Gamma,pF)
%COMPUTES THE DIRICHLET ENERGY (the energy minimized by harmonic maps) of a
%point-to-point map between two shapes
%Src: source shape
%Tar: target shape
%pF: point-to-point map Src -> Tar (a vector of length numvertex(Src) )

%NOTE: the statistical quantitites and the associated terminology do not
%make sense for this energy. Namely:
%mean: really the value of the Dirichlet energy
%standardDeviation: a sort of equivalent of standard deviation, assuming
    %that mean is an expectation value in the quantum-mechanical sense (without normalization)
%worstCase: does not make much sense, either.


    bad_indices = find(~pF);
    pF(bad_indices) = 1;

    D = Src.SHAPE.W .* (Tar_Gamma(pF,pF)/Src.SHAPE.sqrt_area).^2 ;
    
    D(bad_indices,:) = 0;
    D(:,bad_indices) = 0;
    
    
    metrics.mean = abs(0.25*full( sum(sum( D )) ));
    
    metrics.values = [];



end

