function [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_clean_v3(Src_refined, Tar_refined, Steklov_settings)
%PrincipledSteklovBoundaryMatch_clean: Matches the boundary circles by
%matching the normal derivatives of harmonic functions defined by the
%landmarks. This is done through exhaustive search. The normal derivatives
%are to be computed in the ComputeAll function.


toplot = false; %Set to true for plots.


num_shifts = Steklov_settings.num_shifts; % Determines the number of different shifts tried. Includes the trivial shift
% num_it_ICP = Steklov_settings.num_it_ICP;

%% Preliminaries

num_bounds = Steklov_settings.num_landmarks;
num_eigs = Steklov_settings.DS_num_eigs;
nDS = Steklov_settings.num_DS_functions;

ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);

%% Optimize the initial correspondence via a shifting of the orientation

best_shifts = zeros(num_bounds,1); %Stores the answer of the first optimization.

p2pTarSrc = cell(1,num_bounds);
p2pSrcTar = cell(1,num_bounds);

for b = 1:num_bounds % Iterate over each boundary
    
    
    
    Src_nv = length(Src_refined.STEKLOV.boundaries_segment{b});
    Tar_nv = length(Tar_refined.STEKLOV.boundaries_segment{b});
        
%     Src_ND = zeros(Src_nv,num_bounds-1);
%     Tar_ND = zeros(Tar_nv,num_bounds-1);
    
    inds = [1:(b-1) (b+1):num_bounds];
    
    Src_ND = Src_refined.STEKLOV.landmark_harmonics_ND{b}(:, inds);
    Tar_ND = Tar_refined.STEKLOV.landmark_harmonics_ND{b}(:, inds);
    
     
    Tar_line = (1:Tar_nv)/Tar_nv;
    
    bestE = Inf;
    
    for s = 0:(num_shifts - 1)
        
        Src_line = mod( (1:Src_nv)/Src_nv + s/num_shifts , 1);
        
        p2pTarSrc_pre = annquery(Src_line, Tar_line,1);
        p2pSrcTar_pre = annquery(Tar_line, Src_line,1);
        
        Src_on_Tar = Src_ND(p2pTarSrc_pre,:);
        Tar_on_Src = Tar_ND(p2pSrcTar_pre,:);
        
        E = trace( (Tar_ND - Src_on_Tar)' * Tar_refined.STEKLOV.MassBoundary{b} * (Tar_ND - Src_on_Tar) +...
            (Src_ND - Tar_on_Src)' * Src_refined.STEKLOV.MassBoundary{b} * (Src_ND - Tar_on_Src) );       
        
        
        if E < bestE
            
            bestE = E;           
            best_shifts(b) = s; 
%             current_best = s
            
            p2pTarSrc{b} = p2pTarSrc_pre;
            p2pSrcTar{b} = p2pSrcTar_pre;
            
        end


    end
    
   
    
    
end


 
%% Extract the Fmap on harmonic functions from the boundary p2p map

FmapSrcTar = Principled_FmapFromBoundary_P2P(Src_refined, Tar_refined, p2pTarSrc, Steklov_settings);
FmapTarSrc = Principled_FmapFromBoundary_P2P(Tar_refined, Src_refined, p2pSrcTar, Steklov_settings);





if toplot
    figure('Units','normalized','Position',[0.25 0.25 0.5 0.5])

    subplot(2,3,1)
    imagesc(FmapSrcTar); colorbar;
    subplot(2,3,2)
    imagesc(FmapTarSrc); colorbar;
    subplot(2,3,3)
    imagesc(FmapTarSrc*FmapSrcTar); colorbar;
    
    for bc_to_plot = 1:num_bounds
        subplot(2,3,4)
        plot(p2pTarSrc{bc_to_plot},'.'); title(sprintf('p2p Tar->Src, It: %d',0)); hold on
        subplot(2,3,5)
        plot(p2pSrcTar{bc_to_plot},'.'); title(sprintf('p2p Src->Tar, It: %d',0)); hold on
        subplot(2,3,6)
        plot(p2pTarSrc{bc_to_plot}(p2pSrcTar{bc_to_plot}),'.'); title(sprintf('p2p Tar->Tar inv check, It: %d',0)); hold on
    end
    
    
    
end






%% Convert the Fmap to a full p2p map (on segment (non-landmark refined vertices) )

[Src_FTarSrc, Tar_FTarSrc] = build_Principled_EmbeddingsForS2_to_S1_p2p(Src_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Src_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Src_refined.STEKLOV.W_segment,...
                                                                                 Tar_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Tar_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Tar_refined.STEKLOV.W_segment,...
                                                                                 FmapSrcTar, FmapTarSrc, OrthoWeight, ProperWeight, BijectWeight, Steklov_settings);       
                                                                              
[Tar_FSrcTar, Src_FSrcTar] = build_Principled_EmbeddingsForS2_to_S1_p2p(Tar_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Tar_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Tar_refined.STEKLOV.W_segment,...
                                                                                 Src_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Src_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Src_refined.STEKLOV.W_segment,...
                                                                                 FmapTarSrc, FmapSrcTar, OrthoWeight, ProperWeight, BijectWeight, Steklov_settings); 
                   
fullp2pTarSrc = SeparateInteriorBoundaryNearestNeighbor_segment(Src_refined, Src_FTarSrc, Tar_refined, Tar_FTarSrc);
fullp2pSrcTar = SeparateInteriorBoundaryNearestNeighbor_segment(Tar_refined, Tar_FSrcTar, Src_refined, Src_FSrcTar);


if toplot
    
    figure('Units','normalized','Position',[0.25 0.25 0.5 0.5])
    
    subplot(1,2,1)
    plot(fullp2pTarSrc(fullp2pSrcTar), '.')
    subplot(1,2,2)
    plot(fullp2pSrcTar(fullp2pTarSrc), '.')
    
    
    
end











end