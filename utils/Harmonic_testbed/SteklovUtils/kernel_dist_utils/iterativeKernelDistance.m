function [D,dD] = iterativeKernelDistance(E1, E2, sigma, power)
%ITERATIVEKERNELDISTANCE computes the (square) gaussian kernel distance and its derivative for subsequent iterations.
%E1: embedding of shape 1
%E2: embedding of shape 2 --> This embedding is deformed during optimization
%sigma: "standard deviation" of gaussian kernel (kernel is not normalized as probability)
%power: the vector of the powers of the embedding E2 -> variable over which
%optimizations take place


power = power(:)'; % line vector of powers.
E2d = 0.5 * E2 + 0.5 * E2.^(power + 1); % E2 ---> 0.5 * E2 (1 + E2^power) 


K11 = gaussianKernel(E1, E1, sigma);
K12 = gaussianKernel(E1, E2d, sigma);
K22 = gaussianKernel(E2d, E2d, sigma);

D = sum(K11,[1 2]) + sum(K22,[1 2]) - 2 * sum(K12, [1 2]);


%% Gradient

E2der = E2.^power;

dD = zeros(1, length(power) );

for i = 1:length(power)
    
    dK22 = -1/(sigma^2) * power(i) * K22 .* ( E2d(:,i) - E2d(:,i)' ) .* ( E2der(:,i) - E2der(:,i)' );
    
    dK12 = 1/(sigma^2) * power(i) * K12 .* ( E1(:,i) - E2d(:,i)' ) .* ( zeros(size(E1(:,i))) - E2der(:,i)');
    
    
    dD(i) = sum(dK22,[1 2]) - 2 * sum(dK12,[1 2]); 
    
    
    
end

dD = dD(:);


end

