function reorder_lengths = computeLengthAlongBoundary(boundary_reorder, Shape)
% Produces an arclength parametrization of the boundary.

%boundary_reorder: output of boundaryReorder
%Shape: structure with field .surface.VERT

num_boundaries = length(boundary_reorder);

for b = 1:num_boundaries

    reorder_lengths{b} = zeros(length(boundary_reorder{b}), 1);

    for i = 2:length(reorder_lengths{b})

        reorder_lengths{b}(i) = reorder_lengths{b}(i-1) + sqrt(sum( ( Shape.surface.VERT(boundary_reorder{b}(i),:) - Shape.surface.VERT(boundary_reorder{b}(i-1),:)  ).^2 ,2));


    end

end





end