function corr_out = reinsertBoundaryCorrespondence_Principled(p2p_S1S2, S1_landmarks, S2_landmarks, S2_segment)
% Reinserts the known landmark correspondence into a p2p map on segments.
%p2p_S1S2: p2p map from interior S1 to interior S2
%S1_landmarks: list of landmark vertices of S1 in the original indexing
%S2_landmarks: same for S2

corr_out = zeros(1,length(p2p_S1S2) + length(S1_landmarks));

non_lmks = 1:(length(p2p_S1S2) + length(S1_landmarks));
non_lmks(S1_landmarks) = [];


corr_out(non_lmks) = S2_segment(p2p_S1S2); %takes the different orderings into account.
corr_out(S1_landmarks) = S2_landmarks;


end