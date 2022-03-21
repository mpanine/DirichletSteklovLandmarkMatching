function p2p = SeparateInteriorBoundaryNearestNeighbor_segment(S1_refined, S1_F, S2_refined, S2_F)
%function p2p = SeparateInteriorBoundaryNearestNeighbor_segment
%A nearest neighbor search on functions that explicitly maps the 
%the landmark boundaries to the corresponding landmark boundaries
%and the rest to the rest. Works on the "segment" portion of the refined
%shapes.

%The produced p2p map goes S2 -> S1


%% DEBUG

% S1_real = isreal(S1_F)
% S2_real = isreal(S2_F)
% 
% S1_nan = sum(sum(isnan(S1_F)))
% S2_nan = sum(sum(isnan(S2_F)))

%% Preliminaries - Find the 'rest' of the shapes. Could be done elsewhere for efficiency.

S1_rest = 1:size(S1_refined.STEKLOV.FullBasis_segment,1);
S1_rest( vertcat(S1_refined.STEKLOV.boundaries_segment{:}) ) = [];

S2_rest = 1:size(S2_refined.STEKLOV.FullBasis_segment,1);
S2_rest( vertcat(S2_refined.STEKLOV.boundaries_segment{:}) ) = [];


% S1_rest = 1:S1_refined.SHAPE.nv;
% S1_rest([S1_landmarks; S1_refined.STEKLOV.boundary_list]) = [];
% 
% S2_rest = 1:S2_refined.SHAPE.nv;
% S2_rest([S2_landmarks; S2_refined.STEKLOV.boundary_list]) = [];

%% Computing the p2p map

p2p = zeros(1, size(S2_refined.STEKLOV.FullBasis_segment,1) );

p2p(S2_rest) = S1_rest(...
                        annquery(S1_F(S1_rest, :)',...
                        S2_F(S2_rest, :)', 1)...
                       ); %On the rest

for i = 1:length(S1_refined.STEKLOV.boundaries) % On the boundaries
    
   p2p(S2_refined.STEKLOV.boundaries_segment{i}) = S1_refined.STEKLOV.boundaries_segment{i}(...
                                                        annquery(S1_F(S1_refined.STEKLOV.boundaries_segment{i}, :)',...
                                                        S2_F(S2_refined.STEKLOV.boundaries_segment{i}, :)', 1)....
                                                    );  
    
end
                            
% p2p(S2_landmarks) = S1_landmarks; %On the landmarks; exact correspondence.                           
                          




end