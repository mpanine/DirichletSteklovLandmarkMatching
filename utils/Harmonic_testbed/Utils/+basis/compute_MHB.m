function [BASIS] = compute_MHB(S, MHB_settings)
%Compute the eigen vectors of the cotangent Laplacian matrix

assert(isa(MHB_settings, 'basis.MHB_Settings'), ...
    sprintf('MHB_settings has to be of type basis.MHB_Settings. Got %s instead.', ...
    class(MHB_settings)));

numEigs = MHB_settings.numBasisFun;

if nargin < 2, numEigs = 200; end

try
    [evecs, evals] = eigs(S.W, S.A, numEigs, 1e-5);
catch
    % In case of trouble make the laplacian definite
    [evecs, evals] = eigs(S.W - 1e-8*speye(S.nv), S.A, numEigs, 'sm');
end
evals = diag(evals);

[evals, order] = sort(abs(evals),'ascend');
evecs = evecs(:,order);

BASIS.basis = evecs;
BASIS.basis_inverse = BASIS.basis'*S.A;
BASIS.eigenvalues = evals;
end