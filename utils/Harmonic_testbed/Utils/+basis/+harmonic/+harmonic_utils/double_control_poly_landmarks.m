function double_control_poly_landmarks(obj,event_obj,S1,S2,X1,X2,B1,B2)
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
        
        
        funObj = @(v) controlLandmarkDistance(Basis1, Basis2, landmarks1, landmarks2, NumControlLandmarks, v);
        funProj = @(F) F;
        
        %----------------------------------------------> not implemented from here
        %down
        
        
        
        powers = ones( length(landmarks2) - NumControlLandmarks , 1);
                
        options.verbose = 1;
        options.maxIter = 500;
        options.optTol = 1e-15; % 1e-5 default.
        options.progTol = 1e-15; % 1e-9 default.
        options.numDiff = 1; %Acitvates finite differences.
        
        
        [InitialDistance, InitGrad] = controlLandmarkDistance(Basis1, Basis2, landmarks1, landmarks2, NumControlLandmarks, powers);
        InitialDistance
        InitGrad;
        
        powers = minConf_PQN( funObj, powers, funProj, options);
        powers = powers';
        
        [OptimizedDistance, OptGrad] = controlLandmarkDistance(Basis1, Basis2, landmarks1, landmarks2, NumControlLandmarks, powers);
        OptimizedDistance
        OptGrad;
        
        
        powers_back = ones( length(landmarks2)  - NumControlLandmarks  , 1);
        funObj_back = @(v) controlLandmarkDistance(Basis2.^powers, Basis1, landmarks2, landmarks1, NumControlLandmarks, v);
        
        
        powers_back = minConf_PQN( funObj_back, powers_back, funProj, options);
        powers_back = powers_back';
        
        
        
        
        
        
    pF_test = annquery(Basis1', Basis2', 1);
    pF_test_2 = annquery(Basis2', Basis1', 1);
    
    pF_opt = annquery(Basis1', (Basis2.^powers)', 1);
    pF_opt_2 = annquery( (Basis2.^powers)', Basis1', 1);
    
    pF_opt_back = annquery( (Basis1.^powers_back)', (Basis2.^powers)', 1);
    pF_opt_back_2 = annquery( (Basis2.^powers)', (Basis1.^powers_back)', 1);
    
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
    title(sprintf('opt p2p map'));
    
    subplot(3,3,8)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1.surface,50));    
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
    trimesh(S1.surface.TRIV, Basis1(:,1), Basis1(:,2), Basis1(:,3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(:,1).^powers(1), Basis2(:,2).^powers(2), Basis2(:,3).^powers(3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Opt Embeddings')
    
    
    
    
    
    
    
    
    
    
    
    subplot(3,3,9)
    trimesh(S1.surface.TRIV, Basis1(:,1).^powers_back(1), Basis1(:,2).^powers_back(2), Basis1(:,3).^powers_back(3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal;  
    hold on;
    trimesh(S2.surface.TRIV, Basis2(:,1).^powers(1), Basis2(:,2).^powers(2), Basis2(:,3).^powers(3), ...
        'FaceColor','interp', ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k'); axis equal; 
    title('Both Opt_2 Embeddings')
    
    
    
    
    
    figure(4)
    
    [~, top_functions_1] = max(Basis1');
    [~, top_functions_2] = max(Basis2');
    
    [~, opt_top_functions_2] = max((Basis2.^powers)');
    
    [~, back_opt_top_functions_1] = max((Basis1.^powers_back)');
    
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
    trimesh(S2.surface.TRIV, S2.surface.X, S2.surface.Y, S2.surface.Z, ...
        opt_top_functions_2,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Opt Regions S2')
    
    subplot(1,4,4)
    trimesh(S1.surface.TRIV, S1.surface.X, S1.surface.Y, S1.surface.Z, ...
        back_opt_top_functions_1,...
        'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;    
    title('Opt_2 Regions S1')
    
    
    
    
    figure(5)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_test, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_test_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('p2p map'));
    
        
    figure(6)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('opt p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('opt p2p map'));
    
    
    figure(7)
    subplot(2,1,1)    
    visualize_map_lines(S2.surface, S1.surface, ...
        pF_opt_back, ...
        euclidean_fps(S2.surface,50));    
    title(sprintf('opt back p2p map'));
    
    subplot(2,1,2)    
    visualize_map_lines(S1.surface, S2.surface, ...
        pF_opt_back_2, ...
        euclidean_fps(S1.surface,50));    
    title(sprintf('opt back p2p map'));


    
end