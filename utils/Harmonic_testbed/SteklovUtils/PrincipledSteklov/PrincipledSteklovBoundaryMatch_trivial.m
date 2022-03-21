function [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_trivial(Src_refined, Tar_refined, Steklov_settings)
%PrincipledSteklovBoundaryMatch_trivial: Matches the boundary circles in a
%trivial manner: circles go to circles and that is about all there is to
%it. The overall orientation of the circle is given by the outwards normal
%of the shape. No further optimization is done.


toplot = false; %Set to true for plots.


%% Preliminaries

num_bounds = Steklov_settings.num_landmarks;
nDS = Steklov_settings.num_DS_functions;

ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);


%% Build the trivial circle to circle maps

p2pTarSrc = cell(1,num_bounds);
p2pSrcTar = cell(1,num_bounds);


for b = 1:num_bounds % Trivial circle maps on all boundaries.
    
    Src_nv = size(Src_refined.STEKLOV.DS_evecs_bound{b},1);
    Tar_nv = size(Tar_refined.STEKLOV.DS_evecs_bound{b},1);
    
    Tar_line = (1:Tar_nv)/Tar_nv;
    Src_line = (1:Src_nv)/Src_nv;
        
    p2pTarSrc{b} = annquery(Src_line, Tar_line,1);
    p2pSrcTar{b} = annquery(Tar_line, Src_line,1);
        
  
end

%% Convert the circle to circle maps to a Fmap on Dirichlet-Steklov basis

FmapSrcTar = Principled_FmapFromBoundary_P2P(Src_refined, Tar_refined, p2pTarSrc, Steklov_settings);
FmapTarSrc = Principled_FmapFromBoundary_P2P(Tar_refined, Src_refined, p2pSrcTar, Steklov_settings);

if toplot
    figure
    subplot(1,3,1)
    imagesc(FmapSrcTar); colorbar;
    subplot(1,3,2)
    imagesc(FmapTarSrc); colorbar;
    subplot(1,3,3)
    imagesc(FmapSrcTar*FmapSrcTar); colorbar;
end

%% Convert the Fmaps to a full p2p map (on segment)

[Src_FTarSrc, Tar_FTarSrc] = build_Principled_EmbeddingsForS2_to_S1_p2p(Src_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Src_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Src_refined.STEKLOV.W_segment,...
                                                                                 Tar_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Tar_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Tar_refined.STEKLOV.W_segment,...
                                                                                 FmapSrcTar, FmapTarSrc, OrthoWeight, ProperWeight, BijectWeight, Steklov_settings);       
                                                                              
[Tar_FSrcTar, Src_FSrcTar] = build_Principled_EmbeddingsForS2_to_S1_p2p(Tar_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Tar_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Tar_refined.STEKLOV.W_segment,...
                                                                                 Src_refined.STEKLOV.FullBasis_segment(:,1:nDS) , Src_refined.STEKLOV.FullBasisProduct(1:nDS,1:nDS), Src_refined.STEKLOV.W_segment,...
                                                                                 FmapTarSrc, FmapSrcTar, OrthoWeight, ProperWeight, BijectWeight, Steklov_settings); 
                   
fullp2pTarSrc = annquery(Src_FTarSrc' , Tar_FTarSrc', 1);
fullp2pSrcTar = annquery(Tar_FSrcTar' , Src_FSrcTar', 1);  


end