function [D,dD] = controlLandmarkDistance(E1, E2, landmarks1, landmarks2, num_control, power)
%CO0NTROLLANDMARKDISTANCE computes the (square) distance between functions on control landmarks and its derivative for subsequent iterations.
%E1: embedding of shape 1
%E2: embedding of shape 2 --> This embedding is deformed during optimization
%landmarks1: landmarks on E1
%landmarks2: landmarks on E2
%num_control: number of control landmarks (last few landmarks)
%power: the vector of the powers of the embedding E2 -> variable over which
%optimizations take place


power = power(:)'; % line vector of powers.

values1 = E1(landmarks1(end + 1 - num_control:end), : ); %Values at the control landmarks
values2 = E2(landmarks2(end + 1 - num_control:end), : );

uu = values2.^power;
u = values1 - uu;

D = norm( u )^2;


%% Gradient

pre_dD = -2 * u .* log(values2) .* uu;
dD = sum(pre_dD,1);
dD = dD(:);


end
