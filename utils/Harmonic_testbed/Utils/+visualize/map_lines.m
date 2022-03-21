function map_lines(S1,S2,map,samples)

    X1 = S1.surface.VERT(:,1);
    Y1 = S1.surface.VERT(:,2);
    Z1 = S1.surface.VERT(:,3);
    
    X2 = S2.surface.VERT(:,1);
    Y2 = S2.surface.VERT(:,2);
    Z2 = S2.surface.VERT(:,3);

    g1 = normalize_function(0.1,0.99,Y2);
    g2 = normalize_function(0.1,0.99,Z2);
    g3 = normalize_function(0.1,0.99,X2);

    f1 = g1(map);
    f2 = g2(map);
    f3 = g3(map);

    % plot semi-transparent meshes
    trimesh(S1.surface.TRIV, X1, Y1, Z1, ...
        'FaceVertexCData', [f1 f2 f3], 'FaceColor','interp', ...
        'FaceAlpha', 0.9, 'EdgeColor', 'none'); axis equal;
    hold on;
    
    xdiam = 3/2*(max(X2)-min(X2));
    trimesh(S2.surface.TRIV, X2+xdiam, Y2, Z2, ...
        'FaceVertexCData', [g1 g2 g3], 'FaceColor','interp', ...
        'FaceAlpha', 0.9, 'EdgeColor', 'none'); axis equal;
    
    if (~isempty(samples))
        target_samples = map(samples);
        
        Xstart = X1(samples)'; Xend = X2(target_samples)';
        Ystart = Y1(samples)'; Yend = Y2(target_samples)';
        Zstart = Z1(samples)'; Zend = Z2(target_samples)';
        
        Xend = Xend+xdiam;
        Colors = [f1 f2 f3];
        ColorSet = Colors(samples,:);
        set(gca, 'ColorOrder', ColorSet);
        plot3([Xstart; Xend], [Ystart; Yend], [Zstart; Zend]);
    end    
end