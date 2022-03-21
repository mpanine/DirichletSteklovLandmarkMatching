function [ BASIS ] = compute_LDB( S, LDB_settings )
%COMPUTE_LDB computes the diffusion function basis, given input parameters.
    assert(isa(LDB_settings, 'basis.LDB_Settings'), ...
    sprintf('LDB_settings has to be of type basis.LDB_Settings. Got %s instead.', ...
    class(LDB_settings)));

    if ~isfield(S, 'landmarkIndices')
        throw(MException('COMPUTE_LDBERROR:potentialNotDefined',...
                         sprintf('landmarkIndices needs to be defined on the shape %s to allow the computation of the LDB.',...
                                    S.name)));
    end
    diffusionBasis = zeros(S.nv, LDB_settings.numBasisFun);
    nf = LDB_settings.parameters.numFunPerMu;
    for i=1:length(LDB_settings.parameters.mu)
        mu = LDB_settings.parameters.mu(i);
        diffusionBasis(:, (i-1)*nf+1:i*nf) = diffusion_basis_noPrint(S.surface.VERT,...
                                            mu, S.landmarkIndices(1:nf),...
                                            S.W, S.A);
    end;
    
    diffusionBasis = full(diffusionBasis);
    
    % Orthonormalizing the NON-orthogonal diffusion basis using Gram-Schmidt
    [orth_diffusionBasis] = weighted_orthonormal_basis_noPrint(diffusionBasis, S.A);
    
    orth_diffusionBasis_inverse = orth_diffusionBasis'*S.A;
    
    BASIS.basis = orth_diffusionBasis;
    BASIS.basis_inverse = orth_diffusionBasis_inverse;
end

