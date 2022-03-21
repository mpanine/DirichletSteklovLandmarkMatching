        load 0001.null.0.mat;
        
        [W, A] = mshlp_matrix(shape);

        % Eigendecomposition
        Am = sparse([1:length(A)], [1:length(A)], A);
        [evecs evals] = eigs(W, Am, 200, -1e-5);
        evals = diag(evals);

        count = 1;
        for t = [256 512 1024 4096],
            LL = exp(-t*abs(evals));
            K  = evecs*diag(LL(:))*evecs';
            D  = sqrt(repmat(diag(K),[1 size(K,2)]) + repmat(diag(K)',[size(K,2) 1]) - 2*K);
            subplot(1,4,count); count = count+1;
            trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, D(1,:)); axis image; axis off; shading interp;
        end
        
        