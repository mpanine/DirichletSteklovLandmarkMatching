function [S1_opt, harmonic_basis_1_opt, current_radf] = centralFEoptimization(S1, S2, harmonic_basis_1, harmonic_basis_2, landmarks1, landmarks2, NumControlLandmarks, maxradii_1)
%CENTRALFEOPTIMIZATION optimizes the radii of the landmarks in order to
%match the values of the basis functions on control landmarks
%S1: shape being optimized
%S2: The other shape
%harmonic_basis_1: Harmonic basis of the unrefined S1
%harmonic_basis_2: harmonic basis of the unrefiend S2
%landmarks1: landmarks on S1, incuding NumControlLandmarks control landmarks at the end
%landmarks2: landmarks on S2, incuding NumControlLandmarks control landmarks at the end
%NumControlLandmarks: number of control landmarks used


%Internal Settings
minradf = 0.1; %minima and maxima for 
maxradf = 0.5;

current_radf = (maxradf - minradf)/2 * ones(size(maxradii_1));


% num_circles = 2; % Internal setting

control1 = landmarks1(end-NumControlLandmarks+1 : end);
control2 = landmarks2(end-NumControlLandmarks+1 : end);



S1_opt = S1;
harmonic_basis_1_opt = harmonic_basis_1;

it = 0
for l = 1:length(landmarks1)-NumControlLandmarks
    
    stepsign = 1;
    
    while true 
    
        current_dist = norm( harmonic_basis_1_opt(control1, l) - harmonic_basis_2(control2, l) );
%         dist = Inf

        radf = current_radf;
        radf(l) = radf(l) + stepsign*0.05; % First we modify the radius according to stepsign

        if radf(l) > maxradf
           radf(l) = maxradf; 
        end

        S1_opt_t = insertCentralFElist(S1_opt, landmarks1(1:end-NumControlLandmarks), radf.*maxradii_1);
        harmonic_basis_1_opt_t = ComputeCentralHarmonicBasis(S1_opt_t, landmarks1(1:end-NumControlLandmarks));
        
        if norm(imag( harmonic_basis_1_opt_t )) > 0
           break; 
        end
        
        dist = norm( harmonic_basis_1_opt_t(control1, l) - harmonic_basis_2(control2, l) )

        if dist < current_dist %Did we improve the situation by changing the radius?

            current_radf = radf;
            S1_opt = S1_opt_t;
            harmonic_basis_1_opt = harmonic_basis_1_opt_t;
            
            

        else %try changing the sign of the step
                
                stepsign = -stepsign;
                
                radf = current_radf;
                radf(l) = radf(l) + stepsign*0.05;

                if radf(l) < minradf
                   radf(l) = minradf; 
                end

                S1_opt_t = insertCentralFElist(S1_opt, landmarks1(1:end-NumControlLandmarks), radf.*maxradii_1);
                harmonic_basis_1_opt_t = ComputeCentralHarmonicBasis(S1_opt_t, landmarks1(1:end-NumControlLandmarks));
                
                if norm(imag( harmonic_basis_1_opt_t )) > 0
                     break; 
                end

                dist = norm( harmonic_basis_1_opt_t(control1, l) - harmonic_basis_2(control2, l) )    
            
                if dist < current_dist

                    current_radf = radf;
                    S1_opt = S1_opt_t;
                    harmonic_basis_1_opt = harmonic_basis_1_opt_t;
                    
                else %Altering the radius does not help. Stop the iteration.
                    
                    break;
                    
                end
                
                
        end
    
    it= it +1
    end
    

end




end