function [S1_opt, harmonic_basis_1_opt] = refinementOptimization(S1, S2, harmonic_basis_1, harmonic_basis_2, landmarks1, landmarks2, num_control_landmarks)
%REFINEMENTOPTIMIZATION refines the meshes around the landmarks to make the
%harmonic functions on the control landmarks as close to each other as
%possible.
%S1: shape being optimized
%S2: The other shape
%harmonic_basis_1: Harmonic basis of the unrefined S1
%harmonic_basis_2: harmonic basis of the unrefiend S2
%landmarks1: landmarks on S1, incuding num_control_landmarks control landmarks at the end
%landmarks2: landmarks on S2, incuding num_control_landmarks control landmarks at the end
%num_control_landmarks: number of control landmarks used


num_circles = 2; % Internal setting

control1 = landmarks1(end-num_control_landmarks+1 : end);
control2 = landmarks2(end-num_control_landmarks+1 : end);



S1_opt = S1;
harmonic_basis_1_opt = harmonic_basis_1;

it = 0;
for l = 1:length(landmarks1)-num_control_landmarks
    
    current_dist = norm( harmonic_basis_1(control1, l) - harmonic_basis_2(control2, l) );
    dist = Inf;
    
    while true % The optimization continues until no further improvement can be made
    
        if true %sum( harmonic_basis_1(control1, l) < harmonic_basis_2(control2, l) ) ~= num_control_landmarks %the function on shape 1 is not smaller on all landmarks -> optimization can make sense

            S1.surface = refineVertexNeighborhood(S1.surface, landmarks1(l), num_circles);
            S1 = recomputeShape(S1);

            harmonic_basis_1 = ComputeHarmonicBasis(S1, landmarks1(1:end-num_control_landmarks) );
            dist = norm( harmonic_basis_1(control1, l) - harmonic_basis_2(control2, l) );

            if dist < current_dist

                S1_opt = S1;
                harmonic_basis_1_opt = harmonic_basis_1;
                current_dist = dist;
                
            else
                break;
            end
            
        else

            break;
            
        end
        
        it = it + 1

    end

end




end