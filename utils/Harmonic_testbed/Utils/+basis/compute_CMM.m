function [ BASIS ] = compute_CMM( S, CMM_settings )
%COMPUTE_CMM computes the compressed manifold mode basis, given input 
% parameters.
    assert(isa(CMM_settings, 'basis.CMM_Settings'), ...
        sprintf('CMM_settings has to be of type basis.CMM_Settings. Got %s instead.', ...
        class(CMM_settings)));
    
    if CMM_settings.parameters.accelerated
        [basis, ~, ~] = acmm(S.W, S.D, CMM_settings.parameters.mu,...
                                       CMM_settings.numBasisFun);
        basis_inverse = basis'*S.D;
    else
        basis = compressed_manifold_mode(S.surface.VERT,S.surface.TRIV,-S.W, S.D,...
            CMM_settings.numBasisFun, CMM_settings.parameters.mu);
        basis_inverse = basis'*S.D;
    end
    
    BASIS.basis = basis;
    BASIS.basis_inverse = basis_inverse;
end