function [ Energy_F, dEnergy_F, Energy_v, dEnergy_v, x0_v, F_init, est_rank ] = compute_partial( Src, Tar, energy_settings )
    %COMPUTE_PARTIAL Summary of this function goes here
    %   Detailed explanation goes here
    
    numBasisFun_Src = energy_settings.basis_settings_Src.numBasisFun;
    numBasisFun_Tar = energy_settings.basis_settings_Tar.numBasisFun;

    basis_Src = Src.BASIS.basis(:,1:numBasisFun_Src);
%     basis_Tar = Tar.BASIS.basis(:,1:numBasisFun_Tar);
    basis_inverse_Src = Src.BASIS.basis_inverse(1:numBasisFun_Src,:);
    basis_inverse_Tar = Tar.BASIS.basis_inverse(1:numBasisFun_Tar,:);
    eigenValues_Src = Src.SHAPE.laplacianBasis.eigenvalues(1:numBasisFun_Src);
    eigenValues_Tar = Tar.SHAPE.laplacianBasis.eigenvalues(1:numBasisFun_Tar);
    
    mu1 = energy_settings.weight_slantedDiagonalMask;
    mu2 = energy_settings.weight_subOrthogonality;
    mu3 = energy_settings.weight_area;
    mu4 = energy_settings.weight_regularity;
    
    opt.tv_sigma = energy_settings.tv_sigma;
    opt.tv_mean = energy_settings.tv_mean;

    % COMPUTE W MASK
    est_rank = sum(eigenValues_Tar-max(eigenValues_Src)<0);
    W = zeros(numBasisFun_Src);
    for i=1:numBasisFun_Src
        for j=1:numBasisFun_Src
            slope = est_rank/numBasisFun_Src;
            direction = [1 slope];
            direction = direction./norm(direction);
            W(i,j) = exp(-0.03*sqrt(i.^2 + j.^2))*norm(cross([direction 0], [i,j, 0]-[1 1 0]));        
        end
    end
    
    est_rank = sum(Tar.SHAPE.laplacianBasis.eigenvalues...
                            (1:numBasisFun_Tar)...
                            -...
                       max(Src.SHAPE.laplacianBasis.eigenvalues...
                                (1:numBasisFun_Src))...
                                < 0);
                            
        fprintf('Estimated rank: %d\n', est_rank);

    d = zeros(1,numBasisFun_Src);
    d(1:est_rank) = 1;

    D = repmat(d,numBasisFun_Src,1);
    
    F_init = (max(max(W))-W)./max(max(W));
    
    A = basis_inverse_Tar*Tar.DESCRIPTORS.fct;
    B = @(v) basis_inverse_Src*diag(sparse(v))*Src.DESCRIPTORS.fct;

    Energy_F = @(v, F) (...
            sum(sum((F*A-B(v)).^2).^0.5) + ...
            mu1 * norm(F.*W,'fro')^2 + ...
            mu2 * (norm(F'*F,'fro')^2 - sum(diag(F'*F).^2) + ...
                    sum((diag(F'*F) - d').^2) ));

    dEnergy_F = @(v, F) (...
                        norm_21_gradient(F,A,B(v)) + ...
                        mu1 * 2 * F.*W.*W + ...
                        mu2 * 4 * (F*F'*F - F.*D ));
                    
    A_v = @(F) F*basis_inverse_Tar*Tar.DESCRIPTORS.fct;
    B_v = basis_inverse_Src;
    
    x0_v = @(F) basis_Src*(F*(basis_inverse_Tar*ones(Tar.SHAPE.nv,1)));
    
    area_Tar = full( sum(diag(Tar.SHAPE.D)) );
    areas_Src = full(diag(Src.SHAPE.D));
       
    vfunc = mumford_shah(Src.SHAPE.surface.VERT,...
                         Src.SHAPE.surface.TRIV, Src.SHAPE.D);
    Energy_v = @(v, F) vfunc.cost(A_v(F), B_v, Src.DESCRIPTORS.fct, v, ones(size(v)), area_Tar, areas_Src, mu3, mu4, opt);
    dEnergy_v = @(v, F) vfunc.grad(A_v(F), B_v, Src.DESCRIPTORS.fct, v, ones(size(v)), area_Tar, areas_Src, mu3, mu4, opt);
    
end
