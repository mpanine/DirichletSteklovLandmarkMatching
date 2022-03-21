function [D,dD] = ICPowDistance(E1, E2, power, p2p, numskip)
%ICPowDistance computes the (square) distance between functions on images of a p2p map and its derivative for subsequent iterations.
%E1: embedding of shape 1
%E2: embedding of shape 2 --> This embedding is deformed during optimization
%power: powers applied to E2: variables over which 
%p2p: point-to-point map going E2 -> E1
%numskip: Allows to skip some number of points to speed up the code. Set to 1 to have no skip.


power = power(:)'; % line vector of powers.

% values1 = E1(1:numskip:length(p2p), :); %Values at the images of the p2p map
% values2 = E2(p2p(1:numskip:end) , : );

values1 = E1(p2p(1:numskip:end), :); %Values at the images of the p2p map
values2 = E2(1:numskip:length(p2p) , : );

% size(values1)
% size(values2)

% size(power)

uu = values2.^power;
u = values1 - uu;

D = norm( u )^2;


%% Gradient

pre_dD = -2 * u .* log(values2 + 0.0000001) .* uu;
dD = sum(pre_dD,1);
dD = dD(:);


end
