function [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_clean_v2(Src_refined, Tar_refined, Steklov_settings)
%PrincipledSteklovBoundaryMatch_clean_v2: Matches the boundary circles by
%minimizing the orthonomality energy landmark by landmark.
%Ends by extracting the full p2p map between "segments". (non-landmark vertices of the refined shape)

toplot = false; %Set to true for plots.


num_shifts = Steklov_settings.num_shifts; % Determines the number of different shifts tried. Includes the trivial shift

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
          
    inds = [1:(b-1) (b+1):num_bounds];
    
    Tar_line = (1:Tar_nv)/Tar_nv;
    
    bestE = Inf;
    
    if toplot
        figure
    end
    
    for s = 0:(num_shifts - 1)
        
        Src_line = mod( (1:Src_nv)/Src_nv + s/num_shifts , 1);
        
        p2pTarSrc_pre = annquery(Src_line, Tar_line,1);
        p2pSrcTar_pre = annquery(Tar_line, Src_line,1);
        
        Src_on_Tar = Src_refined.STEKLOV.FullBasis( Src_refined.STEKLOV.boundaries{b}(p2pTarSrc_pre), ( (b-1)*num_eigs+1 ):b*num_eigs );
        Tar_on_Src = Tar_refined.STEKLOV.FullBasis( Tar_refined.STEKLOV.boundaries{b}(p2pSrcTar_pre), ( (b-1)*num_eigs+1 ):b*num_eigs );
        
        % Now we compute the normal derivatives of those - necessary for
        % the orthogonality energy. Uses projection onto the basis of the
        % other shape
        
             
        
        Src_on_Tar_ND = Tar_refined.STEKLOV.FullBasis(Tar_refined.STEKLOV.boundaries{b}, ( (b-1)*num_eigs+1 ):b*num_eigs ) * ...
                        Tar_refined.STEKLOV.DS_evals(:,:,b).^(1.5) *  Tar_refined.STEKLOV.FullBasis(Tar_refined.STEKLOV.boundaries{b}, ( (b-1)*num_eigs+1 ):b*num_eigs )' * ...
                        Tar_refined.STEKLOV.MassBoundary{b} * Src_on_Tar;
        
        Tar_on_Src_ND = Src_refined.STEKLOV.FullBasis(Src_refined.STEKLOV.boundaries{b}, ( (b-1)*num_eigs+1 ):b*num_eigs ) * ...
                        Src_refined.STEKLOV.DS_evals(:,:,b).^(1.5) * Src_refined.STEKLOV.FullBasis(Src_refined.STEKLOV.boundaries{b}, ( (b-1)*num_eigs+1 ):b*num_eigs )' * ...
                        Src_refined.STEKLOV.MassBoundary{b} * Tar_on_Src;
        
        
        
        
        E = norm( eye(num_eigs) - Src_on_Tar' * Tar_refined.STEKLOV.MassBoundary{b} * Src_on_Tar_ND , 'fro')^2 + ...
            norm( eye(num_eigs) - Tar_on_Src' * Src_refined.STEKLOV.MassBoundary{b} * Tar_on_Src_ND , 'fro')^2 ;
        
        if toplot
            clf
            subplot(1,3,1)
            imagesc(Src_on_Tar' * Tar_refined.STEKLOV.MassBoundary{b} * Src_on_Tar_ND); colorbar;
            subplot(1,3,2)
            imagesc(Tar_on_Src' * Src_refined.STEKLOV.MassBoundary{b} * Tar_on_Src_ND); colorbar;
            subplot(1,3,3)
            imagesc(Src_on_Tar' * Tar_refined.STEKLOV.MassBoundary{b} * Src_on_Tar_ND * Tar_on_Src' * Src_refined.STEKLOV.MassBoundary{b} * Tar_on_Src_ND); colorbar;
            pause(.1)

        end
        
        
   
        
        
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
%     bc_to_plot = 6;
    

    clf
    subplot(2,3,1)
    imagesc(zeros(size(Src_refined.STEKLOV.Steklov_evals,1))); colorbar; title('Fmap Src->Tar, It: 0')
    subplot(2,3,2)
    imagesc(zeros(size(Src_refined.STEKLOV.Steklov_evals,1))); colorbar; title('Fmap Tar->Src, It: 0')
    subplot(2,3,3)
    imagesc(zeros(size(Src_refined.STEKLOV.Steklov_evals,1))); colorbar; title('Inverse Check, It: 0')
    pause(1e-9)
    
    for bc_to_plot = 1:num_bounds
        subplot(2,3,4)
        plot(p2pTarSrc{bc_to_plot},'.'); title(sprintf('p2p Tar->Src, It: %d',0)); hold on
        subplot(2,3,5)
        plot(p2pSrcTar{bc_to_plot},'.'); title(sprintf('p2p Src->Tar, It: %d',0)); hold on
        subplot(2,3,6)
        plot(p2pTarSrc{bc_to_plot}(p2pSrcTar{bc_to_plot}),'.'); title(sprintf('p2p Tar->Tar inv check, It: %d',0)); hold on
    end
    pause(2)

    
    
end


    
FmapSrcTar = Principled_FmapFromBoundary_P2P(Src_refined, Tar_refined, p2pTarSrc, Steklov_settings);
FmapTarSrc = Principled_FmapFromBoundary_P2P(Tar_refined, Src_refined, p2pSrcTar, Steklov_settings);


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