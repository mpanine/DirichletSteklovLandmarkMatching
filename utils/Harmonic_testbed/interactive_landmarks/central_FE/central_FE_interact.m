function central_FE_interact(obj,event_obj,S1,S2,X1,X2,B1,B2)
%     [inisize, largesize, step, reg_weight, initype] = read_params();
%DOUBLE_CONTROL_UTILS 
%     [inisize, largesize, step, reg_weight, initype] = read_params();
        

        NumControlLandmarks = 1; % The number of last landsmarks that are considered "control" landmarks.
            
        global landmarks;
        
        landmarks1 = landmarks{1};
        landmarks2 = landmarks{2};
    
               
        if min(length(landmarks1),length(landmarks2)) < 3
           error('Less than three landmarks.') 
        end
        
        if min(length(landmarks1),length(landmarks2)) < NumControlLandmarks
           error('Too few landmarks.') 
        end
        
        if min(length(landmarks1),length(landmarks2)) - NumControlLandmarks <=  3
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
        
        
        pF_harmonic = annquery(Basis1', Basis2', 1);
        pF_harmonic_2 = annquery(Basis2', Basis1', 1);
        
        %% Unoptimized central FE -- fraction of the smallest inner radius
        
        maxradii_1 = computeLandmarkMaxRadii(S1, landmarks1(1:end-NumControlLandmarks));
        maxradii_2 = computeLandmarkMaxRadii(S2, landmarks2(1:end-NumControlLandmarks));
        
        maxradii_unopt = min(maxradii_1, maxradii_2) * 0.2; % This has to be a min(...). Don't be confused by the terminology.
        
        
%         figure(70)
        S1_central_unopt = insertCentralFElist(S1, landmarks1(1:end-NumControlLandmarks), maxradii_unopt);
        S2_central_unopt = insertCentralFElist(S2, landmarks2(1:end-NumControlLandmarks), maxradii_unopt);
        
        
        Basis1_central_unopt = ComputeCentralHarmonicBasis(S1_central_unopt, landmarks1(1:end-NumControlLandmarks));
        Basis2_central_unopt = ComputeCentralHarmonicBasis(S2_central_unopt, landmarks2(1:end-NumControlLandmarks));
        
        
        
        pF_unopt = annquery(Basis1_central_unopt', Basis2_central_unopt', 1);
        pF_unopt_2 = annquery( Basis2_central_unopt', Basis1_central_unopt', 1);
        
        
        
        
        %% Optimization 
        
        [S1_central_opt, Basis1_central_opt, final_radf_1] = centralFEoptimization(S1, S2, Basis1_central_unopt, Basis2_central_unopt, landmarks1, landmarks2, NumControlLandmarks, maxradii_1);
        
%         [S2_central_opt, Basis2_central_opt, final_radf_2] = centralFEoptimization(S2, S1_central_opt, Basis2_central_unopt, Basis1_central_opt, landmarks2, landmarks1, NumControlLandmarks, maxradii_2);
        
        S2_central_opt = S2_central_unopt;
        Basis2_central_opt = Basis2_central_unopt;
        
%         sum(isnan(Basis1_central_opt), [1 2])
%         sum(isnan(Basis2_central_opt), [1 2])
%         
%         norm( imag(Basis1_central_opt) )
%         norm( imag(Basis2_central_opt) )
        
        pF_opt_back = annquery( Basis1_central_opt', Basis2_central_opt', 1);
        pF_opt_back_2 = annquery( Basis2_central_opt', Basis1_central_opt', 1);
        
        assignin('base','p2p_temp',pF_opt_back)
        
        
        
        
%         radii1_opt = centralFEoptimizer(S1_central_unopt, Basis2_central_unopt, landmarks1, landmarks2, NumControlLandmarks, maxradii_unopt);
%         
%         S1_central_opt = insertCentralFElist(S1_central_unopt, landmarks1(1:end-NumControlLandmarks), radii1_opt);
%         Basis1_central_opt = ComputeCentralHarmonicBasis(S1_central_opt, landmarks1(1:end-NumControlLandmarks));        
%         
%         radii2_opt = maxradii_unopt;% centralFEoptimizer(S2_central_unopt, Basis1_central_opt, landmarks2, landmarks1, NumControlLandmarks, maxradii_unopt);
%         
%         S2_central_opt = insertCentralFElist(S2_central_unopt, landmarks2(1:end-NumControlLandmarks), radii2_opt);
%         Basis2_central_opt = ComputeCentralHarmonicBasis(S2_central_opt, landmarks2(1:end-NumControlLandmarks)); 
%         
%         
%         pF_opt_back = annquery( Basis1_central_opt', Basis2_central_opt', 1);
%         pF_opt_back_2 = annquery( Basis2_central_opt', Basis1_central_opt', 1);
%         
    
    
    
    
    figure(3);
    subplot(3,3,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_harmonic, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('p2p map'));
    
    subplot(3,3,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_harmonic_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('p2p map'));
    
    
    
    subplot(3,3,4)    
    visualize_map_lines(S2.surface, S1_central_unopt.surface, ...
        pF_unopt, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('central p2p map'));
    
    subplot(3,3,5)    
    visualize_map_lines(S1_central_unopt.surface, S2.surface, ...
        pF_unopt_2, ...
        euclidean_fps(S1_central_unopt.surface,50));    
    title(sprintf('central p2p map'));
    
    
    
    subplot(3,3,7)    
    visualize_map_lines(S2_central_unopt.surface, S1_central_unopt.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('opt p2p map'));
    
    subplot(3,3,8)    
    visualize_map_lines(S1_central_unopt.surface, S2_central_unopt.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1_central_unopt.surface,50));    
    title(sprintf('opt p2p map'));
    
    
    
    
    
    
    
    
    
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
    trimesh(S1_central_unopt.surface.TRIV, Basis1_central_unopt(:,1), Basis1_central_unopt(:,2), Basis1_central_unopt(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(:,1), Basis2(:,2), Basis2(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Opt Embeddings')
    
    
    
    
    
    
    
    
    
    
    
    subplot(3,3,9)
    trimesh(S1_central_unopt.surface.TRIV, Basis1_central_unopt(:,1), Basis1_central_unopt(:,2), Basis1_central_unopt(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2_central_unopt.surface.TRIV, Basis2_central_unopt(:,1), Basis2_central_unopt(:,2), Basis2_central_unopt(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Opt_2 Embeddings')
    
    
    
    
    
    figure(4)
    
    [~, top_functions_1] = max(Basis1');
    [~, top_functions_2] = max(Basis2');
    
    [~, opt_top_functions_1] = max(Basis1_central_unopt');
    
    [~, back_opt_top_functions_2] = max(Basis2_central_unopt');
    
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
    trimesh(S1_central_unopt.surface.TRIV, S1_central_unopt.surface.X, S1_central_unopt.surface.Y, S1_central_unopt.surface.Z, ...
        opt_top_functions_1,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Opt Regions S1')
    
    subplot(1,4,4)
    trimesh(S2_central_unopt.surface.TRIV, S2_central_unopt.surface.X, S2_central_unopt.surface.Y, S2_central_unopt.surface.Z, ...
        back_opt_top_functions_2,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Opt_2 Regions S2')
    
    
    
    
    figure(5)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_harmonic, ...
        euclidean_fps(S2.surface,50),'k');    
    title(sprintf('p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_harmonic_2, ...
        euclidean_fps(S1.surface,50),'k');    
    title(sprintf('p2p map'));
    
        
    figure(6)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1_central_unopt.surface, ...
        pF_unopt, ...
        euclidean_fps(S2.surface,50),'k');    
    title(sprintf('central unopt'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1_central_unopt.surface, S2.surface, ...
        pF_unopt_2, ...
        euclidean_fps(S1.surface,50),'k');    
    title(sprintf('central unopt'));
    
    
    figure(7)
    subplot(2,1,1)    
    visualize_map_lines(S2_central_unopt.surface, S1_central_unopt.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2_central_unopt.surface,50),'k');    
    title(sprintf('central opt'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1_central_unopt.surface, S2_central_unopt.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1_central_unopt.surface,50),'k');    
    title(sprintf('central opt'));

    
    
    
    figure(8) %All the embeddings
    sgtitle('All base Embeddings')
    
    
    numfuncs = length(landmarks2) - NumControlLandmarks;
    numplots = nchoosek(numfuncs,3);
    plotper = ceil(sqrt(numplots));
    
    funcombs = nchoosek(1:numfuncs,3);
    
    for p=1:numplots
       
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
    sgtitle('All Opt Embeddings')
    
    
    numfuncs = length(landmarks2) - NumControlLandmarks;
    numplots = nchoosek(numfuncs,3);
    plotper = ceil(sqrt(numplots));
    
    funcombs = nchoosek(1:numfuncs,3);
    
    for p=1:numplots
       
        subplot(plotper,plotper, p)
        
        trimesh(S1_central_unopt.surface.TRIV, Basis1_central_opt(:, funcombs(p,1)), Basis1_central_opt(:, funcombs(p,2)), Basis1_central_opt(:, funcombs(p,3)), ...
            'FaceColor','interp', ...
            'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;
        hold on;
        trimesh(S2_central_unopt.surface.TRIV, Basis2_central_opt(:, funcombs(p,1)), Basis2_central_opt(:, funcombs(p,2)), Basis2_central_opt(:, funcombs(p,3)), ...
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
    visualize_map_lines(S2_central_opt.surface, S1_central_opt.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2_central_opt.surface,50));    
    title(sprintf('p2p map'));
    
    subplot(2,3,4)    
    visualize_map_lines(S1_central_opt.surface, S2_central_opt.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1_central_opt.surface,50));    
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
    trimesh(S1_central_unopt.surface.TRIV, Basis1_central_unopt(:,1).^pow, Basis1_central_unopt(:,2).^pow, Basis1_central_unopt(:,3).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;  
    hold on;
    trimesh(S2_central_unopt.surface.TRIV, Basis2_central_unopt(:,1).^pow, Basis2_central_unopt(:,2).^pow, Basis2_central_unopt(:,3).^pow, ...
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
    trimesh(S1_central_unopt.surface.TRIV, Basis1_central_unopt(:,4).^pow, Basis1_central_unopt(:,5).^pow, Basis1_central_unopt(:,6).^pow, ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'b'); axis equal;  
    hold on;
    trimesh(S2_central_unopt.surface.TRIV, Basis2_central_unopt(:,4).^pow, Basis2_central_unopt(:,5).^pow, Basis2_central_unopt(:,6).^pow, ...
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
    trimesh(S1_central_unopt.surface.TRIV, S1_central_unopt.surface.X, S1_central_unopt.surface.Y, S1_central_unopt.surface.Z, ...
        opt_top_functions_1,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    hold on
    plot3(S1_central_unopt.surface.VERT(landmarks1(1:end-NumControlLandmarks),1),S1_central_unopt.surface.VERT(landmarks1(1:end-NumControlLandmarks),2),S1_central_unopt.surface.VERT(landmarks1(1:end-NumControlLandmarks),3),'.','Color','r','MarkerSize',20)
    plot3(S1_central_unopt.surface.VERT(landmarks1(end-NumControlLandmarks+1:end),1),S1_central_unopt.surface.VERT(landmarks1(end-NumControlLandmarks+1:end),2),S1_central_unopt.surface.VERT(landmarks1(end-NumControlLandmarks+1:end),3),'.','Color','b','MarkerSize',20)
    title('Regions S1')
    
    
    
    subplot(2,3,6)
    trimesh(S2_central_unopt.surface.TRIV, S2_central_unopt.surface.X, S2_central_unopt.surface.Y, S2_central_unopt.surface.Z, ...
        back_opt_top_functions_2,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    hold on
    plot3(S2_central_unopt.surface.VERT(landmarks2(1:end-NumControlLandmarks),1), S2_central_unopt.surface.VERT(landmarks2(1:end-NumControlLandmarks),2), S2_central_unopt.surface.VERT(landmarks2(1:end-NumControlLandmarks),3),'.','Color','r','MarkerSize',20)
    plot3(S2_central_unopt.surface.VERT(landmarks2(end-NumControlLandmarks+1:end),1), S2_central_unopt.surface.VERT(landmarks2(end-NumControlLandmarks+1:end),2), S2_central_unopt.surface.VERT(landmarks2(end-NumControlLandmarks+1:end),3),'.','Color','b','MarkerSize',20)
    title('Regions S2')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end