%
% This code accompanies the paper:
%
% "Partial Functional Correspondence"
% Rodola, Cosmo, Bronstein, Torsello, Cremers
% Computer Graphics Forum 2016
%
% Please cite the paper above if you use this code in your research.
%
% Written by Emanuele Rodola
%
function [D, matches] = icp_refine_partial(BASIS_Src, BASIS_Tar, r, C_init, max_iters)

% NOTE: here we use ann() to computer approx. nearest neighbors.
% You can use your own code (or matlab's knnsearch) for the same task,
% but make sure you get the direction right!

fprintf('Running ICP...\n');

C_init = C_init(:,1:r);

X = BASIS_Tar.basis_inverse(1:r, :);
Y = BASIS_Src.basis_inverse;
tree = ann('init', Y);

matches = ann('search', tree, C_init*X, 1)';

err = sum( sqrt(sum((C_init*X - Y(:,matches)).^2)) );
err = err / (size(X,2)*size(C_init,1));

fprintf('(0) MSE: %.2e\n', err);

if max_iters == 0
    ann('deinit', tree);
    D = C_init;
    return
end

% Start iterations

D_prev = C_init;
err_prev = err;

for i=1:max_iters
    
    [U,~,V] = svd(X * Y(:,matches)');
    D = U * V(:,1:r)';
    D = D';
    
    matches = ann('search', tree, D*X, 1)';

    err = sum( sqrt(sum((D*X - Y(:,matches)).^2)) );
    err = err / (size(X,2)*size(C_init,1));
    
    fprintf('(%d) MSE: %.2e\n', i, err);
    
    if err > err_prev
        fprintf('Local optimum reached.\n');
        D = D_prev;
        break;
    end
    
    if (err_prev - err) < 5e-6
        fprintf('Local optimum reached.\n');
        break;
    end
    
    err_prev = err;
    D_prev = D;

end

ann('deinit', tree);
end
