function p2p = SeparateInteriorExteriorNearestNeighbor(S1_refined, S1_F, S2_refined, S2_F)
%SeparateInteriorExteriorNearestNeighbor: A nearest neighbor search on
%functions that explicitly maps the interior vertices to the interior
%vertices and the non-interior vertices(boundary and landmarks, usually) to
%the non-interior vertices.

%The produced p2p map goes S2 -> S1


%% Preliminaries - Findn the non-interiors of the shapes. Could be done elsewhere for efficiency.

S1_non_Interior = 1:S1_refined.SHAPE.nv;
S1_non_Interior(S1_refined.STEKLOV.Interior) = [];

S2_non_Interior = 1:S2_refined.SHAPE.nv;
S2_non_Interior(S2_refined.STEKLOV.Interior) = [];

%% Computing the p2p map

p2p = zeros(1,S2_refined.SHAPE.nv);

p2p(S2_non_Interior) = annquery(S1_F(S1_non_Interior, :)',...
                                S2_F(S2_non_Interior, :)', 1);

p2p(S2_refined.STEKLOV.Interior) = annquery(S1_F(S1_refined.STEKLOV.Interior, :)',...
                                S2_F(S2_refined.STEKLOV.Interior, :)', 1);                            




end