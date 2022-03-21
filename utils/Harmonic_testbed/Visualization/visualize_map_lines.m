function visualize_map_lines(S1,S2,map,samples,varargin)

    if nargin<5
       edgecolor = 'none';
    else
       edgecolor = varargin{1};
    end

    hold off;
    X1 = S1.X;
    Y1 = S1.Y;
    Z1 = S1.Z;
    
    X2 = S2.X;
    Y2 = S2.Y;
    Z2 = S2.Z;

    g1 = normalize_function(0.1,0.99,Y2);
    g2 = normalize_function(0.1,0.99,Z2);
    g3 = normalize_function(0.1,0.99,X2);

    f1 = g1(map);
    f2 = g2(map);
    f3 = g3(map);

    % plot semi-transparent meshes
    trimesh(S1.TRIV, X1, Y1, Z1, ...
        'FaceVertexCData', [f1 f2 f3], 'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', edgecolor); axis equal;
    hold on;
    
    xdiam = 3/2*(max(X2)-min(X2));
    trimesh(S2.TRIV, X2+xdiam, Y2, Z2, ...
        'FaceVertexCData', [g1 g2 g3], 'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', edgecolor); axis equal;
    
    view([0 90]);
    axis off;
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