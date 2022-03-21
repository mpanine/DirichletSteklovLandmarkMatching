function shape_out = faceAreaRefinement(shape, area_factor, max_iterations)
%FACEAREAREFINEMENT refines a mesh until all faces have an area smaller
%than the average times area_factor. The Laplacian and such are also
%computed at the end.
%shape: structure with field surface
%area_factor: defines the treshold for refinement by area_factor*mean(area)
%max_iterations: maximal number of refinements.

areas = MESH.calc_tri_areas(shape.surface);
% mea = mean(areas);
area_threshold = mean(areas)*area_factor; % The area threshold is fixed with the original areas. It is not updated during refinement.

i = 0;
while i <=max_iterations
    
    areas = MESH.calc_tri_areas(shape.surface);
%     length(areas);
    [ma,ia] = max(areas);
%     ma/mea
    
    
    if ma > area_threshold
        
        i = i + 1;
        
        shape.surface = refineTriangle(shape.surface, ia);
        shape.nf = size(shape.surface.TRIV,1);
        shape.nv = size(shape.surface.X,1);
    
    else
            
        break;
        
    end

    

end

opt = {'IfComputeLB',true, 'numEigs', 0};%, 'IfComputeGeoDist',false,'IfFindEdge',false,'IfFindNeigh',false};
shape_out = recomputeShape(shape,opt{:});



end

