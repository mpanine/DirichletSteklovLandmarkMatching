function [Shape_refined, segment_TRIV, evecs, evals, boundary_edges_old, boundary_old] = ComputeCentralSteklovAll_clean(Shape, landmarks, Stekov_settings)
%Computes the Steklov basis for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark

warning('In this version, the constant eigenfunction is removed from the basis. This is a quick fix that works for CONNECTED SHAPES ONLY.')
warning('In this version, the eigenvalues are normalized wrt to the largest eigenvalue.')

%CURRENTLY THE CONSTANT EIGENFUNCTION IS REMOVED FROM THE SET !!!


%OUTPUT:
%Shape: shape with the new points 
%evecs: Steklov basis with the new points
%evals: eigenvalues

Shape_refined = Shape;

% nv = size(Shape_refined.SHAPE.surface.VERT,1);
% s_basis = zeros(size(S.W,1),length(landmarks));

[Shape_refined.SHAPE, boundaries] = refineCircleAroundCenter(Shape.SHAPE, landmarks, Stekov_settings.landmarks_radii, Stekov_settings.numsubdivisions); %Insert circles at landmarks.

nv = size(Shape_refined.SHAPE.surface.VERT,1); %Just to make sure the size is correctly set.
Shape_refined.SHAPE.nv = nv;

boundary_old = vertcat(boundaries{:});%cell2mat(boundary); %List form of the boundary indices: loses information about separate segments. TODO: output all.
% boundary_old = boundary_old(:);

[WW, SS, segment_TRIV, boundaries_segment, boundary_edges_refined] = splitMeshLandmarkSteklov_clean(Shape_refined, landmarks, boundaries);


%% Normalize mass matrices and compute MassBoundary 

MassBoundary = cell(1,length(landmarks));

for i = 1:length(landmarks) %Cycle over the landmarks. Each landmark gets the Steklov condition in turn.            
  
    SSS = diag( sum(SS(boundaries_segment{i}, boundaries_segment{i}) , 1));
    normSSS = sum(sum(SSS));

    MassBoundary{i} = SSS/normSSS;
    
    SS(boundaries_segment{i}, boundaries_segment{i}) = SS(boundaries_segment{i}, boundaries_segment{i}) / normSSS; %Normalizes each boundary component.
    
end

Shape_refined.STEKLOV.S = operatorBoundaryReinsert(SS, landmarks); %Outputs the normalized Steklov mass matrix.

% total_boundary_length = sum(sum(SS))

    try
        [evecs, evals] = eigs(WW, SS, Stekov_settings.num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        [evecs, evals] = eigs(WW - 1e-8*speye(size(WW)), SS, Stekov_settings.num_eigs, 'sm');
        warning('The Steklov boundary is likely bad.')
    end
    
%     if min(min(evals)) < 0 %Correct a slightly negative lowest eigenvalue.
%        evals(1,1) = 0; 
%     end

    evals(1,:) = []; %REMOVE THE CONSTANT EIGENFUNCTION AND EIGENVALUE.
    evals(:,1) = [];
    
    evecs(:,1) = [];
    
    prods = evecs'*SS*evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %Normalization, just in case (sometimes necessary) (WHY?)
    evecs = evecs * inv_pr;


%% Determine the interior vertices -- Initial interior vertices

Interior = 1:Shape.SHAPE.nv;
Interior(landmarks) = [];


%% Determine the refined interior vertices (not an output)

segment = 1:Shape_refined.SHAPE.nv;
segment(landmarks) = [];

%% Bases restricted to the boundary components.

for i = 1:length(landmarks) % Both bases restricted to the boundary.
    
    Shape_refined.STEKLOV.evecs_bound{i} = evecs(boundaries_segment{i}, : );  
%     size( Shape_refined.STEKLOV.evecs_bound{i})

end




%% Combining the outputs into one neat structure
Shape_refined.STEKLOV.segment_TRIV = segment_TRIV;
% Shape_refined.STEKLOV.evecs = evecs;
Shape_refined.STEKLOV.evecs = extendSegmentFunctionsIntegrals(evecs, Shape_refined.STEKLOV.evecs_bound, MassBoundary, segment, landmarks, Shape_refined.SHAPE.nv);
Shape_refined.STEKLOV.evals = evals/max(max(evals));
Shape_refined.STEKLOV.boundary_edges_refined = boundary_edges_refined;
Shape_refined.STEKLOV.boundaries_segment = boundaries_segment;
Shape_refined.STEKLOV.MassBoundary = MassBoundary;
Shape_refined.STEKLOV.Interior = Interior;
Shape_refined.STEKLOV.boundaries = boundaries;
% Shape_refined.STEKLOV.normalTriangles = normalTriangles;

Shape_refined.STEKLOV.boundary_list = unique(boundary_edges_refined);
% Shape_refined.STEKLOV.harmonic_basis = ComputeRefinedCentralHarmonicBasis(Shape_refined, Shape_refined.STEKLOV.boundaries, landmarks); %Should no longer be necessary.


%% Compute the Dirichlet LB basis and output it

[LB_evecs, LB_evals] = computeDirichletLB(Shape_refined.SHAPE, [landmarks(:); Shape_refined.STEKLOV.boundary_list(:)], 200);
LB_evecs = LB_evecs * diag( LB_evals.^-0.5 );

Shape_refined.STEKLOV.LB_evecs = LB_evecs;
Shape_refined.STEKLOV.LB_evals = LB_evals;



    
    
    


end