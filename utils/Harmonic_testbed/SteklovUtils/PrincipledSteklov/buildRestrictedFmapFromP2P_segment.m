function F = buildRestrictedFmapFromP2P_segment(S1, S2, p2pS2S1, Steklov_settings, zo)
% Builds the Fmap:S1 -> S2 corresponding to the p2p map p2pS2S1.
% The Fmap is restricted to preserve the subspaces (landmark and
% non-harmonic). F is thus block-diagonal.
% zo is the size of the current ZoomOut (starts at least at the number of harmonic basis functions)

num_landmarks = length(Steklov_settings.numsubdivisions);

F = zeros(zo);

for i = 1:num_landmarks
    
    start_ind = (i-1)*Steklov_settings.DS_num_eigs + 1;
    end_ind = i*Steklov_settings.DS_num_eigs;
    
%     F(start_ind:end_ind, start_ind:end_ind) = S2.STEKLOV.FullBasis_segment(:,start_ind:end_ind)' * S2.STEKLOV.W_segment * S1.STEKLOV.FullBasis_segment(p2pS2S1,start_ind:end_ind); % Bulk inner product -- Lower quality?
    
    F(start_ind:end_ind, start_ind:end_ind) = S2.STEKLOV.DS_evals(:,:,i) * S2.STEKLOV.FullBasis_segment(S2.STEKLOV.boundaries_segment{i},start_ind:end_ind)'... Boundary inner product.
                            * S2.STEKLOV.MassBoundary{i} *...
                            S1.STEKLOV.FullBasis_segment(p2pS2S1(S2.STEKLOV.boundaries_segment{i}),start_ind:end_ind);
    
    
    
end

if zo > num_landmarks*Steklov_settings.DS_num_eigs  % There are some LB functions used.
    
    LB_start = num_landmarks*Steklov_settings.DS_num_eigs + 1;
    
    F(LB_start:zo, LB_start:zo) = ...
        S2.STEKLOV.FullBasis_segment(:, LB_start:zo)' * S2.STEKLOV.W_segment * S1.STEKLOV.FullBasis_segment(p2pS2S1, LB_start:zo);   
    

    
    
end






end