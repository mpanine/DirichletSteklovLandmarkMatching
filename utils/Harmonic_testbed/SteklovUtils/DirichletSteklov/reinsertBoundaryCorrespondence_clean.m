function corr_out = reinsertBoundaryCorrespondence_clean(p2p_S1S2, S1_landmarks, S1_Interior,  S2_landmarks, S2_Interior)
% Reinserts the known landmark correspondence.
%p2p_S1S2: p2p map from interior S1 to interior S2
%S1_landmarks: list of landmark vertices of S1 in the original indexing
%S2_landmarks: same for S2
%S2_Interior: list of interior vetices of S2 in the original indexing.

corr_out = zeros(length(S1_Interior) + length(S1_landmarks), 1);

% non_lmks = 1:(length(p2p_S1S2) + length(S1_landmarks));
% non_lmks(S1_landmarks) = [];


corr_out(S1_Interior) = S2_Interior(p2p_S1S2); %takes the different orderings into account.
corr_out(S1_landmarks) = S2_landmarks;


end