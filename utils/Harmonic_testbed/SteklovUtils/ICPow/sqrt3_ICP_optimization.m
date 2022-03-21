function [harmonic_basis_1_opt, p2p_opt] = sqrt3_ICP_optimization(path, S1, S2, harmonic_basis_1, harmonic_basis_2, landmarks1, landmarks2)
%SQRT3_ICP_OPTIMIZATION intrinsically refines the meshes around the landmarks to make the
%harmonic functions as close to each other as possible. The process is ICP-like.
%S1: shape being optimized
%S2: The other shape
%harmonic_basis_1: Harmonic basis of the unrefined S1
%harmonic_basis_2: harmonic basis of the unrefiend S2
%landmarks1: landmarks on S1, incuding num_control_landmarks control landmarks at the end
%landmarks2: landmarks on S2, incuding num_control_landmarks control landmarks at the end


% num_circles = 2; % Internal setting

max_depth = 50; % Internal setting: maximal number of refinements
stepsize = 1;

% control1 = landmarks1(end-num_control_landmarks+1 : end);
% control2 = landmarks2(end-num_control_landmarks+1 : end);



% S1_opt = S1;
harmonic_basis_1_opt = harmonic_basis_1;

it = 0;
refinement_levels = zeros(length(landmarks1),1);


% for l = length(landmarks1)-num_control_landmarks:-1:1
for l = 1:length(landmarks1)
    
%     current_dist = norm( harmonic_basis_1_opt(control1, :) - harmonic_basis_2(control2, :) );
        
       
    while true % The optimization continues until no further improvement can be made
        
        p2p_opt = annquery(harmonic_basis_1_opt', harmonic_basis_2', 1);        
        current_dist = norm( harmonic_basis_1_opt(p2p_opt, :) - harmonic_basis_2(1:length(p2p_opt), :) );
        
        it = it + 1;
    
%         if true %sum( harmonic_basis_1(control1, l) < harmonic_basis_2(control2, l) ) ~= num_control_landmarks %the function on shape 1 is not smaller on all landmarks -> optimization can make sense
            
%             if l == 6
                refinement_levels(l) = refinement_levels(l) + stepsize;
%             end
            
            harmonic_basis_1 = computeIntrinsicBasis_sqrt3( path, S1.nv, landmarks1(1:end), refinement_levels );

            
            p2p_current = annquery(harmonic_basis_1', harmonic_basis_2', 1);
            dist = norm( harmonic_basis_1(p2p_current, :) - harmonic_basis_2(1:length(p2p_current), :) );

            if dist < current_dist

                harmonic_basis_1_opt = harmonic_basis_1;
                p2p_opt = p2p_current;
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