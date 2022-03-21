function [D,dD] = fourierKernelDistance(E1, E2, sigma, coeffs)
%FOURIERKERNELDISTANCE computes the (square) gaussian kernel distance and its derivative.
%E1: embedding of shape 1
%E2: embedding of shape 2 --> This embedding is deformed during optimization
%sigma: "standard deviation" of gaussian kernel (kernel is not normalized as probability)
%coeffs: the Fourier coefficients that deform the embedding E2 -> variable over which
%optimizations take place size = (num_Basis_funcs, highest_fourier_order)

% tic

weight = 1e3;

E2d = fourierDeform(E2, coeffs);


K11 = gaussianKernel(E1, E1, sigma);
K12 = gaussianKernel(E1, E2d, sigma);
K22 = gaussianKernel(E2d, E2d, sigma);

D = sum(K11,[1 2]) + sum(K22,[1 2]) - 2 * sum(K12, [1 2]) + ...
                            weight * norm(coeffs)^2;

% D = norm(coeffs)^2;

% toc
% tic

%% Gradient

% E2der = E2.^(power-1);

dD = zeros(size(coeffs));

for i = 1:size(coeffs,1)
    
    preK22 = -2/(sigma^2) * K22 .* ( E2d(:,i) - E2d(:,i)' );
    preK12 = -2/(sigma^2) * K12 .* ( E1(:,i) - E2d(:,i)' );
         
    for j = 1:size(coeffs,2)
        
        u = sin(pi * j * E2(:,i));
        dK22 = preK22 .* ( u - u' );
        
        dK12 = preK12 .* ( zeros(size(E1(:,i))) - u' );


        dD(i,j) = sum(dK22,[1 2]) - 2 * sum(dK12,[1 2]); 


    end
end

dD = dD + weight * 2 * coeffs;


dD = dD(:);

% dD = 2 * coeffs;
% 
% dD = dD(:);


% toc

end

