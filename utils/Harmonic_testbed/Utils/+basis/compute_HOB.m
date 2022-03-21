function [ BASIS ] = compute_HOB( S, HOB_settings )
%COMPUTE_HOB computes the Hamiltonian operator eigenfunctions basis, given input parameters.
    assert(isa(HOB_settings, 'basis.HOB_Settings'), ...
    sprintf('HOB_settings has to be of type basis.HOB_Settings. Got %s instead.', ...
    class(HOB_settings)));

    if ~isfield(S, 'potential')
        throw(MException('COMPUTE_HOBERROR:potentialNotDefined',...
                         sprintf('The field potential needs to be defined on the shape %s to allow the computation of the HOB.',...
                                    S.name)));
    end
    
    BASIS = struct;
    [BASIS.basis, BASIS.eigenvalues] = hamiltonian_basis(...
                                        S.W, S.A,HOB_settings.numBasisFun,...
                                        S.potential);
                                    
    BASIS.basis_inverse = BASIS.basis'*S.A;
end



