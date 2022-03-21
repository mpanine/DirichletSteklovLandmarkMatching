function [Basis, evals] = computeIntrinsicDirichletLB(W, A, boundary_indices, numEigs)
%computeIntrinsicDirichletLB computes the Laplace-Beltrami basis subject to
%Dirichlet boundary conditions on the boundary indices.
import basis.harmonic.intrinsicTriangulations.*
import basis.harmonic.harmonic_utils.*

%% Infinite potential version


% prePot = 1000000000 * ones(size(boundary_indices));
% Potential = sparse(boundary_indices, boundary_indices, prePot, size(W,1), size(W,2));
% 
% size(W + Potential)
% size(A)
% 
% 
% try chol(W);
%     disp('W is symmetric positive definite.')
% catch ME
%     disp('W is not symmetric positive definite')
% end
% 
% try chol(W + Potential);
%     disp('W + Potential is symmetric positive definite.')
% catch ME
%     disp('W + Potential is not symmetric positive definite')
% end
% 
% 
% try chol(A);
%     disp('A is symmetric positive definite.')
% catch ME
%     disp('A is not symmetric positive definite')
% end
% 
% 
% [Basis, evals] = eigs(W + Potential, A, numEigs, 1e-6);
% 
% 
% realevecs = isreal(Basis)
% realevals = isreal(evals)
% 
% figure
% plot(evals)





%% BACKUP VERSION



            [Wii, ~, Aii] = boundary_poisson_split_intrinsic(W, A, boundary_indices); % The boundary term Wib shouldn't be necessary with Dirichlet BC.


            boundary_values = zeros(length(boundary_indices),1);

            try
                [evecs, evals] = eigs(Wii, Aii, numEigs, 1e-6);
            catch
                % In case of trouble make the laplacian definite
                [evecs, evals] = eigs(Wii - 1e-8*speye(size(Wii,1)), Aii, numEigs, 'sm');
            %     [evecs, evals] = eigs(Wii, Aii, numEigs, 1e-6);
            %     figure
            %     imagesc(Wii)
            %     [evecs, evals] = eigs(Wii + 100*Aii, Aii, numEigs, 'sm');
            end

%             realevecs = isreal(evecs)
%             realevals = isreal(evals)

            evecs = real(evecs);%A quick fix for some issues.
            evals = real(evals);




            evals = diag(evals);

            [evals, order] = sort(abs(evals),'ascend');


            evecs = evecs(:,order);

            Basis = zeros(size(Wii,1) + length(boundary_indices), numEigs);

            for i = 1:numEigs
                Basis(:,i) = BoundaryReinsert(evecs(:,i), boundary_indices, boundary_values );
            end






end