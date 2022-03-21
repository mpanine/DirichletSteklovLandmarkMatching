function Shape_refined = ComputeCentralDirichletSteklovAll_clean(Shape, landmarks, Steklov_settings) 
%ComputeCentralDirichletSteklovBasis_clean: Refines the shape at the landmarks and computes all the relevant function bases.
%This is a cleaned up version used for large scale testing.

%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark
%num_eigs: number of eigenvalues 

%OUTPUT:
%Shape_refined: shape with the new points as well as a new field .STEKLOV
%with the bases, eigenvalues etc




Shape_refined = Shape;

[Shape_refined.SHAPE, boundaries, normalTriangles] = refineCircleAroundCenter(Shape.SHAPE, landmarks, Steklov_settings.landmarks_radii, Steklov_settings.numsubdivisions); %Insert circles at landmarks.
 % boundaries is a cell array of DISCONNECTED COMPONENTS.


nv = size(Shape_refined.SHAPE.surface.VERT,1); %Just to make sure the size is correctly set.
Shape_refined.SHAPE.nv = nv;

[W, S, segment_TRIV, boundaries_segment, boundary_edges_refined] = splitMeshLandmarkDirichletSteklov(Shape_refined, landmarks, boundaries);

Shape_refined.SHAPE.S = operatorBoundaryReinsert(S, landmarks);



evecs = zeros(nv-length(landmarks), Steklov_settings.num_eigs, length(landmarks) );
evals = zeros(Steklov_settings.num_eigs, Steklov_settings.num_eigs, length(landmarks) );


MassBoundary = cell(1,length(landmarks));


for i = 1:length(landmarks) %Cycle over the landmarks. Each landmark gets the Steklov condition in turn.
            
    Dirichlet_boundary = vertcat( boundaries_segment{[1:(i-1) (i+1):length(landmarks)]} ); %All boundary indices, but that of the ith boundary, which is Steklov.
    
    WW = W; %Create the relevant operators.
    SS = S;    
    
    WW(Dirichlet_boundary,:) = [];
    WW(:,Dirichlet_boundary) = [];
    
    SS(Dirichlet_boundary,:) = [];
    SS(:,Dirichlet_boundary) = [];
        
    
    normSS = sum(sum(SS));
    SS = SS/normSS; %Normalizes Steklov boundary length to 1. SUPER IMPORTANT
    
%     MassBoundary{i} = diag(diag( sum(S(boundaries_segment{i}, boundaries_segment{i}) , 1)) )/ normSS; %Lumped mass on boundary component.
    MassBoundary{i} = diag( sum(S(boundaries_segment{i}, boundaries_segment{i}) , 1)) / normSS; %Lumped mass on boundary component.
    
    try
        [pre_evecs, pre_evals] = eigs(WW, SS, Steklov_settings.num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        [pre_evecs, pre_evals] = eigs(WW - 1e-8*speye(size(WW)), SS, Steklov_settings.num_eigs, 'sm');
        warning('The Steklov boundary is likely bad.')
    end

    prods = pre_evecs'*SS*pre_evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %Normalization, just in case (sometimes necessary) (WHY?)
    
    pre_evecs = pre_evecs * inv_pr;
    
    sign_of_first_ef = sign( sum( sign( pre_evecs(:,1) ) ) ); % Fix a positive sign for the first eigenfunction.
    pre_evecs(:,1) = sign_of_first_ef*pre_evecs(:,1);
    
    evecs(:,:,i) = BoundaryReinsertBatch(pre_evecs, Dirichlet_boundary, zeros(size(Dirichlet_boundary)));
    evals(:,:,i) = pre_evals/max(max(pre_evals)); %Normalization by the largest value. TODO: test other normalizations.
  
    
    
end


%% Determine the interior vertices -- Initial interior vertices

Interior = 1:Shape.SHAPE.nv;
Interior(landmarks) = [];


%% Determine the refined interior vertices (not an output)

segment = 1:Shape_refined.SHAPE.nv;
segment(landmarks) = [];

%% Bases restricted to the boundary components.

for i = 1:length(landmarks) % Both bases restricted to the boundary.
    
    Shape_refined.STEKLOV.evecs_bound{i} = evecs(boundaries_segment{i}, : ,i);

end







%% Combining the outputs into one neat structure
Shape_refined.STEKLOV.segment_TRIV = segment_TRIV;
% Shape_refined.STEKLOV.evecs = evecs;
Shape_refined.STEKLOV.evecs = extendSegmentFunctionsIntegrals(evecs, Shape_refined.STEKLOV.evecs_bound, MassBoundary, segment, landmarks, Shape_refined.SHAPE.nv);
Shape_refined.STEKLOV.evals = evals;
Shape_refined.STEKLOV.boundary_edges_refined = boundary_edges_refined;
Shape_refined.STEKLOV.boundaries_new = boundaries_segment;
Shape_refined.STEKLOV.MassBoundary = MassBoundary;
Shape_refined.STEKLOV.Interior = Interior;
Shape_refined.STEKLOV.boundaries = boundaries;
Shape_refined.STEKLOV.normalTriangles = normalTriangles;

Shape_refined.STEKLOV.boundary_list = unique(boundary_edges_refined);
Shape_refined.STEKLOV.harmonic_basis = ComputeRefinedCentralHarmonicBasis(Shape_refined, Shape_refined.STEKLOV.boundaries, landmarks);


%% Compute the Dirichlet LB basis and output it

[LB_evecs, LB_evals] = computeDirichletLB(Shape_refined.SHAPE, [landmarks(:); Shape_refined.STEKLOV.boundary_list(:)], 200);
LB_evecs = LB_evecs * diag( LB_evals.^-0.5 );

Shape_refined.STEKLOV.LB_evecs = LB_evecs;
Shape_refined.STEKLOV.LB_evals = LB_evals;



end