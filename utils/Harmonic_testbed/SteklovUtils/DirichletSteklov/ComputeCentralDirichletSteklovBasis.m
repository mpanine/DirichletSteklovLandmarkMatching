function [Shape_refined, segment_TRIV, evecs, evals, boundary_edges_old, boundaries_new, MassBoundary, Interior, boundaries, normalTriangles, W_segment] = ComputeCentralDirichletSteklovBasis(Shape, landmarks, radii, numsubdivisions, num_eigs)
%Computes the Dirichlet-Steklov basis for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark
%num_eigs: number of eigenvalues 

%OUTPUT:
%Shape: shape with the new points 
%evecs: Dirichlet-Steklov basis with the new points
%evals: eigenvalues
%Interior: List of interior vertices (Only the original ones, no refined vertices included)


% s_basis = zeros(size(S.W,1),length(landmarks));

Shape_refined = Shape;

[Shape_refined.SHAPE, boundaries, normalTriangles] = refineCircleAroundCenter(Shape.SHAPE, landmarks, radii, numsubdivisions); %Insert circles at landmarks.
 % boundaries is a cell array of DISCONNECTED COMPONENTS.

% for i = 1:length(boundaries)
%     
%    iiii  = boundaries{i} 
% end

nv = size(Shape_refined.SHAPE.surface.VERT,1);
Shape_refined.SHAPE.nv = nv;

[W_segment, S, segment_TRIV, boundaries_new, boundary_edges_old] = splitMeshLandmarkDirichletSteklov(Shape_refined, landmarks, boundaries);

Shape_refined.SHAPE.S = operatorBoundaryReinsert(S, landmarks);

% for i = 1:length(boundaries_new)
%     
%    boundaries_new{i} 
% end

evecs = zeros(nv-length(landmarks), num_eigs, length(landmarks) );
evals = zeros(num_eigs, num_eigs, length(landmarks) );


% size(boundaries{1})
% size(boundaries_new{1})



for i = 1:length(landmarks) %Cycle over the landmarks. Each landmark gets the Steklov condition in turn.
    
    Dirichlet_boundary = vertcat( boundaries_new{[1:(i-1) (i+1):length(landmarks)]} ); %All boundary indices, but that of the ith boundary, which is Steklov.
%     size(Dirichlet_boundary)
    
%     component_list = [1:(i-1) (i+1):length(landmarks)];
%     Dirichlet_boundary = exctractBoundaryListFromCellArray(boundaries_new, component_list);
%     size(Dirichlet_boundary)
    
    WW = W_segment; %Create the relevant operators.
    SS = S;
    
    
    WW(Dirichlet_boundary,:) = [];
    WW(:,Dirichlet_boundary) = [];
    
    SS(Dirichlet_boundary,:) = [];
    SS(:,Dirichlet_boundary) = [];
    
    
    
    normSS = sum(sum(SS));
    SS = SS/normSS; %Normalizes Steklov boundary length to 1. SUPER IMPORTANT
    
    MassBoundary{i} = diag( sum(S(boundaries_new{i}, boundaries_new{i}) , 1)) / normSS; %Lumped mass on boundary component.
    
    try
        [pre_evecs, pre_evals] = eigs(WW, SS, num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        [pre_evecs, pre_evals] = eigs(WW - 1e-8*speye(size(WW)), SS, num_eigs, 'sm');
        warning('The Steklov boundary is likely bad.')
    end

    prods = pre_evecs'*SS*pre_evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %Normalization, just in case (sometimes necessary) (WHY?)
    
    pre_evecs = pre_evecs * inv_pr;
    
    sign_of_first_ef = sign( sum( sign( pre_evecs(:,1) ) ) ); % Fix a positive sign for the first eigenfunction.
    pre_evecs(:,1) = sign_of_first_ef*pre_evecs(:,1);
    
    evecs(:,:,i) = BoundaryReinsertBatch(pre_evecs, Dirichlet_boundary, zeros(size(Dirichlet_boundary)));
    evals(:,:,i) = pre_evals;
    

    
end







%% Determine the interior vertices

Interior = 1:Shape.SHAPE.nv;

Interior(landmarks) = [];

















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