function Fmap = Principled_FmapFromBoundary_P2P(Src_refined, Tar_refined, p2p_Tar_Src, Steklov_settings)
%Principled_FmapFromBoundary_P2P extracts functional maps from boundary p2p maps.
%Used in the principled approach.

%p2p: Tar -> Src

% nDS = Steklov_settings.num_DS_functions;

Fmap = zeros(Steklov_settings.num_DS_functions);

for i = 1:Steklov_settings.num_landmarks
    
    start_ind = (i-1)*Steklov_settings.DS_num_eigs + 1;
    end_ind = i*Steklov_settings.DS_num_eigs;

    Fmap(start_ind:end_ind, start_ind:end_ind) = Tar_refined.STEKLOV.DS_evals(:,:,i) *... % Normalization.
                                                 Tar_refined.STEKLOV.FullBasis_segment(Tar_refined.STEKLOV.boundaries_segment{i}, start_ind:end_ind)' * ...
                                                 Tar_refined.STEKLOV.MassBoundary{i} * ...
                                                 Src_refined.STEKLOV.FullBasis_segment(Src_refined.STEKLOV.boundaries_segment{i}(p2p_Tar_Src{i}), start_ind:end_ind);
                                          

end


end