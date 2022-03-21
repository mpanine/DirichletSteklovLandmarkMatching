function [ objectiveFunction, objectiveGradient ] = compute_objectiveFunction( Energies, dEnergies )
%COMPUTE_OBJECTIVEFUNCTION Computes the objective function
%   Inputs:
%       - Energies: cell array of the energies
%       - dEnergies: cell array of the energies derivatives

    objectiveFunction = @(F) sum(cellfun(@(E) E(F), Energies));
    objectiveGradient = @(F) sum(cell2mat(cellfun(@(dE) dE(F), ...
                                 dEnergies, 'UniformOutput', false)), 2);

end

