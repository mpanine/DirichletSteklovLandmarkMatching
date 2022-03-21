function plotCentralHarmonic(Shape, landmarks, funnumber)
%PLOTCENTRALHARMONIC plots the harmonic functions for different internal
%radius sizes


maxradii = computeLandmarkMaxRadii(Shape, landmarks);

thingy = 1:length(landmarks);
maxradii( thingy ~= funnumber) = 0.0000001; % Only one landmark changes radius

Basis0 = ComputeHarmonicBasis(Shape, landmarks );

figure
subplot(3,4,1)
trimesh(Shape.surface.TRIV, Shape.surface.X, Shape.surface.Y, Shape.surface.Z, ...
        Basis0(:,funnumber), 'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
colorbar
title('Harmonic')



radfracs = [0.1:0.1:0.9];


for i = 1:length(radfracs)

    Shape = insertCentralFElist(Shape, landmarks, maxradii*radfracs(i));
    Basis = ComputeCentralHarmonicBasis(Shape, landmarks);
        
    subplot(3,4,i+1)    
    
    trimesh(Shape.surface.TRIV, Shape.surface.X, Shape.surface.Y, Shape.surface.Z, ...
        abs(Basis(:,funnumber) - Basis0(:,funnumber)), 'FaceColor','interp', ...
        'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal;
    colorbar
    
    
    if min(Basis(:,funnumber)) < 0 || max(Basis(:,funnumber)) > 1
        
        title(sprintf('R frac = %0.1f ; MPV = %0.3f', radfracs(i), max(abs(min(Basis(:,funnumber))), abs(max(Basis(:,funnumber)) - 1 ) ) ))
    
    else    
        title(sprintf('R frac = %0.1f', radfracs(i)))
    end




end




end