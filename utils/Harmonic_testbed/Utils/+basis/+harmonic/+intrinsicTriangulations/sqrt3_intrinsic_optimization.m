function harmonic_basis_1_opt = sqrt3_intrinsic_optimization(objpath1, S1, S2, harmonic_basis_1, harmonic_basis_2, landmarks1, landmarks2, num_control_landmarks)
%SQRT#_INTRINSIC_REFINEMENT_OPTIMIZATION intrinsically refines the meshes around the landmarks to make the
%harmonic functions on the control landmarks as close to each other as
%possible.
%S1: shape being optimized
%S2: The other shape
%harmonic_basis_1: Harmonic basis of the unrefined S1
%harmonic_basis_2: harmonic basis of the unrefiend S2
%landmarks1: landmarks on S1, incuding num_control_landmarks control landmarks at the end
%landmarks2: landmarks on S2, incuding num_control_landmarks control landmarks at the end
%num_control_landmarks: number of control landmarks used


% num_circles = 2; % Internal setting

max_depth = 25; % Internal setting: maximal number of refinements
stepsize = 1;

control1 = landmarks1(end-num_control_landmarks+1 : end);
control2 = landmarks2(end-num_control_landmarks+1 : end);



% S1_opt = S1;
harmonic_basis_1_opt = harmonic_basis_1;

it = 0;
refinement_levels = zeros(length(landmarks1) - num_control_landmarks,1);


% for l = length(landmarks1)-num_control_landmarks:-1:1
for l = 1:length(landmarks1)-num_control_landmarks
    
%     current_dist = norm( harmonic_basis_1_opt(control1, :) - harmonic_basis_2(control2, :) );
        
       
    while true % The optimization continues until no further improvement can be made
        
        current_dist = norm( harmonic_basis_1_opt(control1, :) - harmonic_basis_2(control2, :) );
        
        it = it + 1;
    
%         if true %sum( harmonic_basis_1(control1, l) < harmonic_basis_2(control2, l) ) ~= num_control_landmarks %the function on shape 1 is not smaller on all landmarks -> optimization can make sense
            
%             if l == 6
                refinement_levels(l) = refinement_levels(l) + stepsize;
%             end
            
            harmonic_basis_1 = computeIntrinsicBasis_sqrt3( objpath1, S1.nv, landmarks1(1:end-num_control_landmarks), refinement_levels );
            dist = norm( harmonic_basis_1(control1, :) - harmonic_basis_2(control2, :) );

            if dist < current_dist

                harmonic_basis_1_opt = harmonic_basis_1;
                current_dist = dist;
                
                if refinement_levels(l) >= max_depth
                    
                   break;
                    
                end
                
            else
%                 if l == 6
                refinement_levels(l) = refinement_levels(l) - stepsize;
%                 end
                
                break;
            end
            
%         else
% 
%             break;
%             
%         end
        
        

    end

end

% current_dist


end