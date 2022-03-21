function [opt_p2p_S1_to_S2, opt_p2p_S2_to_S1] = isometricZoomOut(S1, S2, Basis1, Basis2, LB1, LB2, p2p_S1_to_S2, p2p_S2_to_S1, start_size, end_size, step_size, landmarks1, landmarks2)

%LB1, LB2: The Laplace-Beltrami basis used. May be either the usual one, or
%the one with Dirichlet boundary conditions.

% maxFailedTries = 50; % A way to stop the iteration earlier. Sets the maximal number of consecutive failed tries.
% FailedTries = 0;


% ifDirichletBasis = true; % Set to true to use the LB basis with Dirichlet boundary conditions.
ifPowerDecay = false; %Set to true to decay the importance of the harmonic part. (This is a bad idea.)


opt_p2p_S1_to_S2 = p2p_S1_to_S2;
opt_p2p_S2_to_S1 = p2p_S2_to_S1;

% figure('position',[200 200 2000 1000])

% best_dist_1 = euclideanSelfDistortion(S1, opt_p2p_S2_to_S1(opt_p2p_S1_to_S2));
% best_dist_2 = euclideanSelfDistortion(S2, opt_p2p_S1_to_S2(opt_p2p_S2_to_S1));


% if ifDirichletBasis
%     
%     [LB1, ~] = computeDirichletLB(S1, landmarks1, end_size);
%     [LB2, ~] = computeDirichletLB(S2, landmarks2, end_size); 
%     
% else
%     
%     LB1 = S1.evecs;
%     LB2 = S2.evecs;
%     
% end

% if ~ifDirichletBasis
%     
%     [LB1, ~] = computeDirichletLB(S1, landmarks1, end_size);
%     [LB2, ~] = computeDirichletLB(S2, landmarks2, end_size); 
%     
% else
%     
%     LB1 = S1.evecs;
%     LB2 = S2.evecs;
%     
% end




power = 0;

for i = start_size:step_size:end_size

%     C12 = LB2(:,1:i) \ LB1(opt_p2p_S2_to_S1,1:i);
%     C21 = LB1(:,1:i) \ LB2(opt_p2p_S1_to_S2,1:i);
    
    C12 = LB2(:,1:i)' * S2.A * LB1(opt_p2p_S2_to_S1,1:i);
    C21 = LB1(:,1:i)' * S1.A * LB2(opt_p2p_S1_to_S2,1:i);
    
    B1 = [Basis1'; LB1(:, 1:i)'];
    B2 = [Basis2'; LB2(:,1:i)'];
    
    if ifPowerDecay
    
        power = power + 1;
        
        BB1 = [( Basis2.^power )'; C21*LB2(:,1:i)'];
        BB2 = [( Basis1.^power )'; C12*LB1(:,1:i)'];
        
        
    else
        
        weight = 1; %sqrt( i/size(Basis2,2) ); % Should compensate for the increase in the size of the LB basis

        BB1 = [weight*Basis2'; C21*LB2(:,1:i)'];
        BB2 = [weight*Basis1'; C12*LB1(:,1:i)'];
    
    end
    
%     BB1 = [weight*Basis2'; C21*(C12*C12')*LB2(:,1:i)'];
%     BB2 = [weight*Basis1'; C12*(C21*C21')*LB1(:,1:i)'];

    opt_p2p_S1_to_S2 = annquery(BB1, B1, 1);
    opt_p2p_S2_to_S1 = annquery(BB2, B2, 1);
    
%     current_dist_1 = euclideanSelfDistortion(S1, pre_opt_p2p_S2_to_S1(pre_opt_p2p_S1_to_S2));
%     current_dist_2 = euclideanSelfDistortion(S2, pre_opt_p2p_S1_to_S2(pre_opt_p2p_S2_to_S1));
    
    
%     clf
%     sgtitle(['Fmap size = ' num2str(i)])
%     subplot(2,2,1)
%     plot(pre_opt_p2p_S2_to_S1(pre_opt_p2p_S1_to_S2),'.')
%     title('S1 to S1')
%     subplot(2,2,2)
%     plot(pre_opt_p2p_S1_to_S2(pre_opt_p2p_S2_to_S1),'.')
%     title('S2 to S2')
%         
%     subplot(2,2,3)
%     toplot = sqrt( sum((S1.surface.VERT(pre_opt_p2p_S2_to_S1(pre_opt_p2p_S1_to_S2),:) - S1.surface.VERT).^2, 2) );
%     plot(toplot, '.')
%     ylim([0 1])    
%     title(['S1 self-error, distortion = ' num2str(current_dist_1) ', best distortion = ' num2str(best_dist_1)])
%     
%     subplot(2,2,4)
%     toplot = sqrt( sum((S2.surface.VERT(opt_p2p_S1_to_S2(opt_p2p_S2_to_S1),:) - S2.surface.VERT).^2, 2) );
%     plot(toplot, '.')
%     ylim([0 1])
%     title(['S2 self-error, distortion = ' num2str(current_dist_2) ', best distortion = ' num2str(best_dist_2)])
    
%     if (current_dist_1 < best_dist_1)&&(current_dist_2 < best_dist_2)
%         opt_p2p_S1_to_S2 = pre_opt_p2p_S1_to_S2;
%         opt_p2p_S2_to_S1 = pre_opt_p2p_S2_to_S1;
%         
%         best_dist_1 = current_dist_1;
%         best_dist_2 = current_dist_2;
%         
%     end
    
%     if current_dist_1 + current_dist_2 < best_dist_1 + best_dist_2
%         opt_p2p_S1_to_S2 = pre_opt_p2p_S1_to_S2;
%         opt_p2p_S2_to_S1 = pre_opt_p2p_S2_to_S1;
%         
%         best_dist_1 = current_dist_1;
%         best_dist_2 = current_dist_2;
%         
%         FailedTries = 0;
%         
%     else
%         
%         FailedTries = FailedTries + 1;
%         
%         if FailedTries >= maxFailedTries
%            break 
%         end
%         
%     end
%     
%     
%     pause(0.0000000000001)
    
end






end