function [Shape, segment_TRIV, evecs, evals, boundary_edges_old, boundary_old] = ComputeCentralSteklovBasis(Shape, landmarks, radii, numsubdivisions, num_eigs)
%Computes the Steklov basis for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark

%OUTPUT:
%Shape: shape with the new points 
%evecs: Steklov basis with the new points
%evals: eigenvalues

nv = size(Shape.SHAPE.surface.VERT,1);
% s_basis = zeros(size(S.W,1),length(landmarks));

[Shape.SHAPE, boundary] = refineCircleAroundCenter(Shape.SHAPE, landmarks, radii, numsubdivisions); %Insert circles at landmarks.

boundary_old = vertcat(boundary{:});%cell2mat(boundary); %List form of the boundary indices: loses information about separate segments. TODO: output all.
% boundary_old = boundary_old(:);

[WW, SS, segment_TRIV, boundary_new, boundary_edges_old] = splitMeshLandmarkSteklov(Shape, landmarks, boundary_old);

% total_boundary_length = sum(sum(SS))

    try
        [evecs, evals] = eigs(WW, SS, num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        [evecs, evals] = eigs(WW - 1e-8*speye(size(WW)), SS, num_eigs, 'sm');
        warning('The Steklov boundary is likely bad.')
    end

    prods = evecs'*SS*evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %Normalization, just in case (sometimes necessary) (WHY?)
    evecs = evecs * inv_pr;




% boundary_long = [];
% boundary_section_length = zeros(size(landmarks));
% 
% for i = 1:length(boundary)
%     boundary_long = [boundary_long; boundary{i}];
%     boundary_section_length(i) = length(boundary{i});
% end
% 
% boundary_section_length = boundary_section_length(:);
% boundary_section_ends = [0; cumsum(boundary_section_length)];
% 
% [W, S, segment_TRIV, boundary_old, boundary_new, boundary_edges_old] = splitSegmentSteklov(S.SHAPE, segment);


% [Wii, Wib, ~] = boundary_poisson_split(S, boundary_long);


% for i = 1:length(landmarks)
%     
%     landmark_values = zeros(length(boundary_long),1);    
%     landmark_values(boundary_section_ends(i) + 1: boundary_section_ends(i+1) ) = 1;
%     
%     pre_h_basis = BoundaryReinsert(mldivide(Wii, -Wib*landmark_values), boundary_long, landmark_values );
%     s_basis(:,i) = pre_h_basis(1:nv); %Restrict to original vertices
%     
%     
% end


end