function [Fmap, p2p] = DirichletSteklovZoomOutV1(Src_evecs_bound, Src_evals, Src_Mass_Boundary, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary)
%SteklovZoomOutV1: First version of the Dirichlet-Steklov ZoomOut

%evecs_bound: eigenvectors on the boundary (cell array over boundary components)
%evals: eigenvalues as a 3d array (diagonal matrices for each boundary component)

% max_size = 20;

Src_evals_mod = Src_evals;
Tar_evals_mod = Tar_evals;


%     for i = 1:length(Src_evecs_bound) % Optional resolvent formulation. DOES NOT WORK RIGHT NOW.
%         Src_evals_mod(:,:,i) = diag( diag(Src_evals(:,:,i) + 10).^-1  );
%         Tar_evals_mod(:,:,i) = diag( diag(Tar_evals(:,:,i) + 10).^-1  );
%     end

%% Guess initial map. Version 1.
% 
% for i = 1:length(Src_evecs_bound) %Guess the initial Fmap
% 
%     Src_guess{i} = sum(Src_Mass_Boundary{i}*Src_evecs_bound{i},1); %Unit function on relevant boundary expanded in basis.
%     Tar_guess{i} = sum(Tar_Mass_Boundary{i}*Tar_evecs_bound{i},1);
%     Fmap{i} = Tar_guess{i}(1,1)' * Src_guess{i}(1,1);
% 
% end
% 
% p2p = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals_mod, Tar_evecs_bound, Tar_evals_mod, Fmap, 1);


%% Guess initial map. Version 2: assume circular-ish correspondence

for i = 1:1:length(Src_evecs_bound)
    
    Src_nv = size(Src_evecs_bound{i},1);
    Tar_nv = size(Tar_evecs_bound{i},1);
    
    Src_line = (1:Src_nv)/Src_nv;
    Tar_line = (1:Tar_nv)/Tar_nv;
    p2p{i} = annquery(Src_line,Tar_line,1);

end

%%

figure

for it = 1:size(Src_evals_mod,2)
    
    Fmap = SteklovDirichletFmapFromP2P(Src_evecs_bound, Tar_evecs_bound, Tar_Mass_Boundary, p2p, it);
    p2p = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals_mod, Tar_evecs_bound, Tar_evals_mod, Fmap, it);    
    
    clf
    imagesc(Fmap{1}); colorbar; title(sprintf('Iteration %d',it))
    pause(1e-9)
    
end






















end