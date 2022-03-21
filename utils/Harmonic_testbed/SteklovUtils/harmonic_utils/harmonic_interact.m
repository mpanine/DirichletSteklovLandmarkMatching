function harmonic_interact(obj,event_obj,S1,S2,X1,X2,B1,B2)
%     [inisize, largesize, step, reg_weight, initype] = read_params();
    
            
        global landmarks;        
        global same_landmarks;
        
        landmarks2 = landmarks{2};
        
        if same_landmarks %Copy the landmarks from shape 1.
            
            landmarks1 = landmarks{2};
            
        else %Use the actual landmarks from shape 2.
            
            landmarks1 = landmarks{1};
           
        end
        
    
        % The code below uses a very simple initialization using landmark
        % preservation with argmin_{C} ||CF - G|| + reg_weight M .* C, 
        % where F and G are delta functions of the landmarks and M is a 
        % matrix of ones everywhere except the diagonal.
        
        if min(length(landmarks1),length(landmarks2)) < 3
           error('Less than three landmarks.') 
        end
        
        if(length(landmarks1) == length(landmarks2))
            
            Basis1 = ComputeHarmonicBasis(S1, landmarks1)';
            Basis2 = ComputeHarmonicBasis(S2, landmarks2)';
            
        else
            warning('uneven number of landmarks %d vs %d',length(landmarks1),...
                length(landmarks2));
        end
    
    pF_test = annquery(Basis1, Basis2, 1);
    pF_test_2 = annquery(Basis2, Basis1, 1);
        
    figure(3);
    subplot(2,3,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_test, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('p2p map'));
    
    subplot(2,3,4)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_test_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('p2p map'));
    
%     subplot(1,3,2)
%     % plot semi-transparent meshes
%     trimesh(S1.surface.TRIV, Basis1(1,:), Basis1(2,:), Basis1(3,:), ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;    
%     title('Embedding S1')
%     
%     subplot(1,3,3)
%     % plot semi-transparent meshes
%     trimesh(S2.surface.TRIV, Basis2(1,:), Basis2(2,:), Basis2(3,:), ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;    
%     title('Embedding S2')
    
%     subplot(2,3,2)
%     
    pow = 1;
%     
%     % plot semi-transparent meshes
%     trimesh(S1.surface.TRIV, Basis1(1,:).^pow, Basis1(2,:).^pow, Basis1(3,:).^pow, ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;    
%     title('Embedding S1')
%     
%     subplot(2,3,3)
%     % plot semi-transparent meshes
%     trimesh(S2.surface.TRIV, Basis2(1,:).^pow, Basis2(2,:).^pow, Basis2(3,:).^pow, ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;    
%     title('Embedding S2')
    
    
    subplot(2,3,2)
    trimesh(S1.surface.TRIV, Basis1(1,:).^pow, Basis1(2,:).^pow, Basis1(3,:).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(1,:).^pow, Basis2(2,:).^pow, Basis2(3,:).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'r'); axis equal; 
    title('Both Embeddings')
    xlabel('$\psi_1$','Interpreter','latex','FontSize',20)
    ylabel('$\psi_2$','Interpreter','latex','FontSize',20)
    zlabel('$\psi_3$','Interpreter','latex','FontSize',20)
    [~, top_functions_1] = max(Basis1);
    [~, top_functions_2] = max(Basis2);
    view([1 1 0.5])
    
    subplot(2,3,3)
    trimesh(S1.surface.TRIV, Basis1(4,:).^pow, Basis1(5,:).^pow, Basis1(6,:).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(4,:).^pow, Basis2(5,:).^pow, Basis2(6,:).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'r'); axis equal; 
    title('Both Embeddings')
    xlabel('$\psi_4$','Interpreter','latex','FontSize',20)
    ylabel('$\psi_5$','Interpreter','latex','FontSize',20)
    zlabel('$\psi_6$','Interpreter','latex','FontSize',20)
%     legend('Shape 1', 'Shape 2')
    [~, top_functions_1] = max(Basis1);
    [~, top_functions_2] = max(Basis2);
    view([1 1 0.5])
    
    
    
    
    subplot(2,3,5)
    trimesh(S1.surface.TRIV, S1.surface.X, S1.surface.Y, S1.surface.Z, ...
        top_functions_1,...
        'FaceColor','flat', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    hold on
    plot3(S1.surface.VERT(landmarks1,1),S1.surface.VERT(landmarks1,2),S1.surface.VERT(landmarks1,3),'.','Color','r','MarkerSize',20)
    title('Regions S1')
    
    
    
    subplot(2,3,6)
    trimesh(S2.surface.TRIV, S2.surface.X, S2.surface.Y, S2.surface.Z, ...
        top_functions_2,...
        'FaceColor','flat', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    hold on
    plot3(S2.surface.VERT(landmarks2,1), S2.surface.VERT(landmarks2,2), S2.surface.VERT(landmarks2,3),'.','Color','r','MarkerSize',20)
    title('Regions S2')
    
    
    
    
%     subplot(1,3,2)
%     % plot semi-transparent meshes
%     trimesh(S1.surface.TRIV, log(Basis1(1,:)), log(Basis1(2,:)), log(Basis1(3,:)), ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;    
%     title('Embedding S1')
%     
%     subplot(1,3,3)
%     % plot semi-transparent meshes
%     trimesh(S2.surface.TRIV, log(Basis2(1,:)), log(Basis2(2,:)), log(Basis2(3,:)), ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;    
%     title('Embedding S2')
    

%     figure(4)
%     
%     trimesh(S1.surface.TRIV, Basis1(1,:).^pow, Basis1(2,:).^pow, Basis1(3,:).^pow, ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
%     hold on;
%     trimesh(S2.surface.TRIV, Basis2(1,:).^pow, Basis2(2,:).^pow, Basis2(3,:).^pow, ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
%     title('Both Embeddings')


    
end