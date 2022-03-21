function [D,dD] = controlLandmarkPolynomialDistance(E1, E2, landmarks1, landmarks2, num_control, coeffs)
%CO0NTROLLANDMARKPOLYNOMIALDISTANCE computes the (square) distance between functions on control landmarks and its derivative for subsequent iterations.
%E1: embedding of shape 1
%E2: embedding of shape 2 --> This embedding is deformed during optimization
%landmarks1: landmarks on E1
%landmarks2: landmarks on E2
%num_control: number of control landmarks (last few landmarks)
%coeffs: the matrix of the polynomial coeffcients of the embedding E2 -> variable over which
%optimizations take place  (SIZE: num_functions * poly_order)



values1 = E1(landmarks1(end + 1 - num_control:end), : ); %Values at the control landmarks
values2 = E2(landmarks2(end + 1 - num_control:end), : );

D = 0;

for c = 1:num_control %Iterate over control landmarks
    
    
    pows = powersOfVector(values2(c,:), size(coeffs,2) );    
    D = D + norm( values1(c,:) - sum(coeffs.^pows, 2) )^2;

    
end


%% Gradient   ---- NOT FINISHED

% dD = zeros(size(coeffs));
% 
% for nu = 1:size(coeffs,1)
%     for mu = 1:size(coeffs,2)
%         
%         for c = 1:num_control
%                 
%             
%             pows = powersOfVector(values2(c,:), size(coeffs,2) );
%             dD(nu,mu) = dD(nu,mu) + ( values1(c,:) - sum(coeffs.^pows, 2) ) * values1(c,)
%             
%         end
%         
%     end
% end
% 
% 
% dD = dD(:);


end
