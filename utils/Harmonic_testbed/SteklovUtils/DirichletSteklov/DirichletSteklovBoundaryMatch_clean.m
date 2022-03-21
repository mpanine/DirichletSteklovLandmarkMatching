function [FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovBoundaryMatch_clean(Src_refined, Tar_refined, Steklov_settings)
%DirichletSteklovBoundaryMatch_clean: Matches the boundary circles by
%matching the normal derivatives of harmonic functions defined by the
%landmarks. Implements Dirichlet-to-Neumann commutativity, orthogonality and
%bijectivity

%THIS VERSION IMPLEMENTS AN EXHAUSTIVE SEARCH FOR THE OPTIMAL INITIAL
%GUESS BASED ON NORMAL DERIVATIVE MATCHING. THEN IT PERFORMS A FEW ICP-LIKE ITERATIONS TO SMOOTH EVERYTHING OUT.


toplot = false; %Set to true for plots.


num_shifts = Steklov_settings.num_shifts; % Determines the number of different shifts tried. Includes the trivial shift
num_it_ICP = Steklov_settings.num_it_ICP;

%% Preliminaries

num_bounds = length(Src_refined.STEKLOV.evecs_bound);
num_eigs = size(Src_refined.STEKLOV.evals,2);

%% Optimize the initial correspondence via a shifting of the orientation

best_shifts = zeros(num_bounds,1); %Stores the answer of the first optimization.

p2pTarSrc = cell(1,num_bounds);
p2pSrcTar = cell(1,num_bounds);

for b = 1:num_bounds % Iterate over each boundary
    
    
    
    Src_nv = size(Src_refined.STEKLOV.evecs_bound{b},1);
    Tar_nv = size(Tar_refined.STEKLOV.evecs_bound{b},1);
        
    Src_ND = zeros(Src_nv,num_bounds-1);
    Tar_ND = zeros(Tar_nv,num_bounds-1);
    
    inds = [1:(b-1) (b+1):num_bounds];
    
    for i = 1:(num_bounds-1)
         
        Src_ND(:,i) = normalDerivativeAtLandmarkBoundaryBoundaryTriangles(Src_refined.SHAPE.surface, Src_refined.STEKLOV.normalTriangles{b}, Src_refined.STEKLOV.harmonic_basis(:,inds(i)), Src_refined.STEKLOV.boundaries{b}, Src_refined.STEKLOV.MassBoundary{b} );
        Tar_ND(:,i) = normalDerivativeAtLandmarkBoundaryBoundaryTriangles(Tar_refined.SHAPE.surface, Tar_refined.STEKLOV.normalTriangles{b}, Tar_refined.STEKLOV.harmonic_basis(:,inds(i)), Tar_refined.STEKLOV.boundaries{b}, Tar_refined.STEKLOV.MassBoundary{b} );
        
    end
    
    
    
    Tar_line = (1:Tar_nv)/Tar_nv;
    
    bestE = Inf;
    
    for s = 0:(num_shifts - 1)
        
        Src_line = mod( (1:Src_nv)/Src_nv + s/num_shifts , 1);
        
        p2pTarSrc_pre = annquery(Src_line, Tar_line,1);
        p2pSrcTar_pre = annquery(Tar_line, Src_line,1);
        
        Src_on_Tar = Src_ND(p2pTarSrc_pre,:);
        Tar_on_Src = Tar_ND(p2pSrcTar_pre,:);
        
        E = (Tar_ND - Src_on_Tar)' * Tar_refined.STEKLOV.MassBoundary{b} * (Tar_ND - Src_on_Tar) +...
            (Src_ND - Tar_on_Src)' * Src_refined.STEKLOV.MassBoundary{b} * (Src_ND - Tar_on_Src);       
        
        
        if E < bestE
            
            bestE = E;            
            best_shifts(b) = s; 
%             current_best = s
            
            p2pTarSrc{b} = p2pTarSrc_pre;
            p2pSrcTar{b} = p2pSrcTar_pre;
            
        end


    end
    
   
    
    
end


 






if toplot
    figure('Units','normalized','Position',[0.25 0.25 0.5 0.5])
    bc_to_plot = 2;
    

    clf
    subplot(2,3,1)
    imagesc(zeros(size(Src_refined.STEKLOV.evals,1))); colorbar; title('Fmap Src->Tar, It: 0')
    subplot(2,3,2)
    imagesc(zeros(size(Src_refined.STEKLOV.evals,1))); colorbar; title('Fmap Tar->Src, It: 0')
    subplot(2,3,3)
    imagesc(zeros(size(Src_refined.STEKLOV.evals,1))); colorbar; title('Inverse Check, It: 0')
    pause(1e-9)
    
    subplot(2,3,4)
    plot(p2pTarSrc{bc_to_plot},'.'); title(sprintf('p2p Tar->Src, It: %d',0))
    subplot(2,3,5)
    plot(p2pSrcTar{bc_to_plot},'.'); title(sprintf('p2p Src->Tar, It: %d',0))
    subplot(2,3,6)
    plot(p2pTarSrc{bc_to_plot}(p2pSrcTar{bc_to_plot}),'.'); title(sprintf('p2p Tar->Tar inv check, It: %d',0))
    pause(2)
    
    

    
    
end












for it = 1:num_it_ICP
    
    FmapSrcTar = SteklovDirichletFmapFromP2P(Src_refined.STEKLOV.evecs_bound, Tar_refined.STEKLOV.evecs_bound, Tar_refined.STEKLOV.MassBoundary, p2pTarSrc, num_eigs);
    FmapTarSrc = SteklovDirichletFmapFromP2P(Tar_refined.STEKLOV.evecs_bound, Src_refined.STEKLOV.evecs_bound, Src_refined.STEKLOV.MassBoundary, p2pSrcTar, num_eigs);
%     p2pTarSrc = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals_mod, Tar_evecs_bound, Tar_evals_mod, FmapSrcTar, it);    
    
    [p2pTarSrc, p2pSrcTar] = DirichletSteklovFmapToP2P_clean(Src_refined.STEKLOV.evecs_bound, Src_refined.STEKLOV.evals,...
                                                             Tar_refined.STEKLOV.evecs_bound, Tar_refined.STEKLOV.evals,...
                                                             FmapSrcTar, FmapTarSrc, Steklov_settings);
    
    
    
    if toplot
        clf
        subplot(2,3,1)
        imagesc(FmapSrcTar{bc_to_plot}); colorbar; title(sprintf('Fmap Src->Tar, It: %d',it))
        subplot(2,3,2)
        imagesc(FmapTarSrc{bc_to_plot}); colorbar; title(sprintf('Fmap Tar->Src, It: %d',it))
        subplot(2,3,3)
        imagesc(FmapSrcTar{bc_to_plot}*FmapTarSrc{bc_to_plot}); colorbar; title(sprintf('Inverse Check, It: %d',it))
        pause(1e-9)
        
        subplot(2,3,4)
        plot(p2pTarSrc{bc_to_plot},'.'); title(sprintf('p2p Tar->Src, It: %d',it))
        subplot(2,3,5)
        plot(p2pSrcTar{bc_to_plot},'.'); title(sprintf('p2p Src->Tar, It: %d',it))
        subplot(2,3,6)
        plot(p2pTarSrc{bc_to_plot}(p2pSrcTar{bc_to_plot}),'.'); title(sprintf('p2p Tar->Tar inv check, It: %d',it))
        pause(.5)
        
        
    end
    
end






















end