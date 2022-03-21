function double_intrinsic_sqrt3(obj,event_obj,S1,S2,X1,X2,B1,B2)
%DOUBLE_CONTROL_UTILS 
%     [inisize, largesize, step, reg_weight, initype] = read_params();
        

        NumControlLandmarks = 1; % The number of last landsmarks that are considered "control" landmarks.
            
        global landmarks;        
        global same_landmarks;
        
        landmarks2 = landmarks{2};
        
        if same_landmarks %Copy the landmarks from shape 1.
            
            landmarks1 = landmarks{2};
            
        else %Use the actual landmarks from shape 2.
            
            landmarks1 = landmarks{1};
           
        end
    
               
        if min(length(landmarks1),length(landmarks2)) < 3
           error('Less than three landmarks.') 
        end
        
        if min(length(landmarks1),length(landmarks2)) < NumControlLandmarks
           error('Too few landmarks.') 
        end
        
        if min(length(landmarks1),length(landmarks2)) - NumControlLandmarks < 3 %<=  3
           error('Too few landmarks.') 
        end
        
        if(length(landmarks1) == length(landmarks2))
            
            Basis1 = ComputeHarmonicBasis(S1, landmarks1(1:end-NumControlLandmarks) );
            Basis2 = ComputeHarmonicBasis(S2, landmarks2(1:end-NumControlLandmarks) );
            
%             size(Basis1)
%             size(Basis2)
            
        else
            warning('uneven number of landmarks %d vs %d',length(landmarks1),...
                length(landmarks2));
        end
        
        
        pF_test = annquery(Basis1', Basis2', 1);
        pF_test_2 = annquery(Basis2', Basis1', 1);
        
        
%         [S1, Basis1_opt] = refinementOptimization(S1, S2, Basis1, Basis2, landmarks1, landmarks2, NumControlLandmarks);
        Basis1_opt = sqrt3_intrinsic_optimization(1, S1, S2, Basis1, Basis2, landmarks1, landmarks2, NumControlLandmarks);
        
        
        pF_opt = annquery(Basis1_opt', Basis2', 1);
        pF_opt_2 = annquery( Basis2', Basis1_opt', 1);
        
        
%         [S2, Basis2_opt] = refinementOptimization(S2, S1, Basis2, Basis1_opt, landmarks2, landmarks1, NumControlLandmarks);
        Basis2_opt = sqrt3_intrinsic_optimization(2, S2, S1, Basis2, Basis1_opt, landmarks2, landmarks1, NumControlLandmarks);
        
        pF_opt_back = annquery( Basis1_opt', Basis2_opt', 1);
        pF_opt_back_2 = annquery( Basis2_opt', Basis1_opt', 1);
        
    
    
    
    
    figure(3);
    subplot(3,3,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_test, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('p2p map'));
    
    subplot(3,3,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_test_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('p2p map'));
    
    
    
    subplot(3,3,4)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('opt p2p map'));
    
    subplot(3,3,5)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('opt p2p map'));
    
    
    
    subplot(3,3,7)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('both opt map'));
    
    subplot(3,3,8)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('both opt map'));
    
    
    
    
    
    
    
    
    
    subplot(3,3,3)
    trimesh(S1.surface.TRIV, Basis1(:,1), Basis1(:,2), Basis1(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(:,1), Basis2(:,2), Basis2(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Embeddings')
    
    
    subplot(3,3,6)
    trimesh(S1.surface.TRIV, Basis1_opt(:,1), Basis1_opt(:,2), Basis1_opt(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(:,1), Basis2(:,2), Basis2(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Opt Embeddings')
    
    
    
    
    
    
    
    
    
    
    
    subplot(3,3,9)
    trimesh(S1.surface.TRIV, Basis1_opt(:,1), Basis1_opt(:,2), Basis1_opt(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2_opt(:,1), Basis2_opt(:,2), Basis2_opt(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Opt_2 Embeddings')
    
    
    
    
    
    figure(4)
    
    [~, top_functions_1] = max(Basis1');
    [~, top_functions_2] = max(Basis2');
    
    [~, opt_top_functions_1] = max(Basis1_opt');
    
    [~, back_opt_top_functions_2] = max(Basis2_opt');
    
    subplot(1,4,1)
    trimesh(S1.surface.TRIV, S1.surface.X, S1.surface.Y, S1.surface.Z, ...
        top_functions_1,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Regions S1')
    
    
    
    subplot(1,4,2)
    trimesh(S2.surface.TRIV, S2.surface.X, S2.surface.Y, S2.surface.Z, ...
        top_functions_2,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Regions S2')
    
    subplot(1,4,3)
    trimesh(S1.surface.TRIV, S1.surface.X, S1.surface.Y, S1.surface.Z, ...
        opt_top_functions_1,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Opt Regions S1')
    
    subplot(1,4,4)
    trimesh(S2.surface.TRIV, S2.surface.X, S2.surface.Y, S2.surface.Z, ...
        back_opt_top_functions_2,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Opt_2 Regions S2')
    
    
    
    
    figure(5)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_test, ...
        euclidean_fps(S2.surface,50),'k');    
    title(sprintf('p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_test_2, ...
        euclidean_fps(S1.surface,50),'k');    
    title(sprintf('p2p map'));
    
        
    figure(6)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt, ...
        euclidean_fps(S2.surface,50),'k');    
    title(sprintf('opt p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_2, ...
        euclidean_fps(S1.surface,50),'k');    
    title(sprintf('opt p2p map'));
    
    
    figure(7)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2.surface,50),'k');    
    title(sprintf('opt back p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1.surface,50),'k');    
    title(sprintf('opt back p2p map'));

    
    
    
    figure(8) %All the embeddings
    sgtitle('All base Embeddings')
    
    
    numfuncs = length(landmarks2) - NumControlLandmarks;
    numplots = nchoosek(numfuncs,3);
    plotper = ceil(sqrt(numplots));
    
    funcombs = nchoosek(1:numfuncs,3);
    
    for p=1:min(numplots,30)
       
        subplot(plotper,plotper, p)
        
        trimesh(S1.surface.TRIV, Basis1(:, funcombs(p,1)), Basis1(:, funcombs(p,2)), Basis1(:, funcombs(p,3)), ...
            'FaceColor','interp', ...
            'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;
        hold on;
        trimesh(S2.surface.TRIV, Basis2(:, funcombs(p,1)), Basis2(:, funcombs(p,2)), Basis2(:, funcombs(p,3)), ...
            'FaceColor','interp', ...
            'FaceAlpha', 0.6, 'EdgeColor', 'r'); axis equal;
%         title('Both Opt_2 Embeddings')
        xlabel(num2str(funcombs(p,1)))
        ylabel(num2str(funcombs(p,2)))
        zlabel(num2str(funcombs(p,3)))
              
        
    end
    
    
    
    
    figure(9) %All the embeddings
    sgtitle('All Opt_2 Embeddings')
    
    
    numfuncs = length(landmarks2) - NumControlLandmarks;
    numplots = nchoosek(numfuncs,3);
    plotper = ceil(sqrt(numplots));
    
    funcombs = nchoosek(1:numfuncs,3);
    
    for p=1:min(numplots,30)
       
        subplot(plotper,plotper, p)
        
        trimesh(S1.surface.TRIV, Basis1_opt(:, funcombs(p,1)), Basis1_opt(:, funcombs(p,2)), Basis1_opt(:, funcombs(p,3)), ...
            'FaceColor','interp', ...
            'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;
        hold on;
        trimesh(S2.surface.TRIV, Basis2_opt(:, funcombs(p,1)), Basis2_opt(:, funcombs(p,2)), Basis2_opt(:, funcombs(p,3)), ...
            'FaceColor','interp', ...
            'FaceAlpha', 0.6, 'EdgeColor', 'r'); axis equal;
%         title('Both Opt_2 Embeddings')
        xlabel(num2str(funcombs(p,1)))
        ylabel(num2str(funcombs(p,2)))
        zlabel(num2str(funcombs(p,3)))
              
        
    end
        
        
    
    %% Plot for report
    
    figure(10);
    subplot(2,3,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('p2p map'));
    
    subplot(2,3,4)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_back_2, ...
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
    trimesh(S1.surface.TRIV, Basis1_opt(:,1).^pow, Basis1_opt(:,2).^pow, Basis1_opt(:,3).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2_opt(:,1).^pow, Basis2_opt(:,2).^pow, Basis2_opt(:,3).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'r'); axis equal; 
    title('Both Embeddings')
    xlabel('$\psi_1$','Interpreter','latex','FontSize',20)
    ylabel('$\psi_2$','Interpreter','latex','FontSize',20)
    zlabel('$\psi_3$','Interpreter','latex','FontSize',20)
%     [~, top_functions_1] = max(Basis1_opt);
%     [~, top_functions_2] = max(Basis2_opt);
    view([1 1 0.5])
    
    subplot(2,3,3)
    trimesh(S1.surface.TRIV, Basis1_opt(:,4).^pow, Basis1_opt(:,5).^pow, Basis1_opt(:,6).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2_opt(:,4).^pow, Basis2_opt(:,5).^pow, Basis2_opt(:,6).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'r'); axis equal; 
    title('Both Embeddings')
    xlabel('$\psi_4$','Interpreter','latex','FontSize',20)
    ylabel('$\psi_5$','Interpreter','latex','FontSize',20)
    zlabel('$\psi_6$','Interpreter','latex','FontSize',20)
%     legend('Shape 1', 'Shape 2')
%     [~, top_functions_1] = max(Basis1_opt);
%     [~, top_functions_2] = max(Basis2_opt);
    view([1 1 0.5])
    
    
    
    
    subplot(2,3,5)
    trimesh(S1.surface.TRIV, S1.surface.X, S1.surface.Y, S1.surface.Z, ...
        opt_top_functions_1,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    hold on
    plot3(S1.surface.VERT(landmarks1(1:end-NumControlLandmarks),1),S1.surface.VERT(landmarks1(1:end-NumControlLandmarks),2),S1.surface.VERT(landmarks1(1:end-NumControlLandmarks),3),'.','Color','r','MarkerSize',20)
    plot3(S1.surface.VERT(landmarks1(end-NumControlLandmarks+1:end),1),S1.surface.VERT(landmarks1(end-NumControlLandmarks+1:end),2),S1.surface.VERT(landmarks1(end-NumControlLandmarks+1:end),3),'.','Color','b','MarkerSize',20)
    title('Regions S1')
    
    
    
    subplot(2,3,6)
    trimesh(S2.surface.TRIV, S2.surface.X, S2.surface.Y, S2.surface.Z, ...
        back_opt_top_functions_2,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    hold on
    plot3(S2.surface.VERT(landmarks2(1:end-NumControlLandmarks),1), S2.surface.VERT(landmarks2(1:end-NumControlLandmarks),2), S2.surface.VERT(landmarks2(1:end-NumControlLandmarks),3),'.','Color','r','MarkerSize',20)
    plot3(S2.surface.VERT(landmarks2(end-NumControlLandmarks+1:end),1), S2.surface.VERT(landmarks2(end-NumControlLandmarks+1:end),2), S2.surface.VERT(landmarks2(end-NumControlLandmarks+1:end),3),'.','Color','b','MarkerSize',20)
    title('Regions S2')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
end