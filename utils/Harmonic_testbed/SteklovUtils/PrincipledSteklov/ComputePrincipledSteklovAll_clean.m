function [Shape_refined, segment_TRIV, Steklov_evecs, Steklov_evals, boundary_edges_old, boundary_old] = ComputePrincipledSteklovAll_clean(Shape, landmarks, Steklov_settings)
%Computes the both the Steklov and Dirichlet-Steklov bases for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark

warning('In this version, the constant eigenfunction is removed from the basis. This is a quick fix that works for CONNECTED SHAPES ONLY.')
warning('In this version, the eigenvalues are normalized wrt to the largest eigenvalue.')

%OUTPUT:
%Shape: shape with the new points 
%evecs: Steklov basis with the new points
%evals: eigenvalues

Shape_refined = Shape;

% nv = size(Shape_refined.SHAPE.surface.VERT,1);
% s_basis = zeros(size(S.W,1),length(landmarks));

[Shape_refined.SHAPE, boundaries] = refineCircleAroundCenter(Shape.SHAPE, landmarks, Steklov_settings.landmarks_radii, Steklov_settings.numsubdivisions); %Insert circles at landmarks.

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

%% Compute the Steklov basis, except the constant eigenfunction

    try
        [Steklov_evecs, Steklov_evals] = eigs(WW, SS, Steklov_settings.Steklov_num_eigs + 1, 1e-5); % The +1 takes care of the removed constant eigenfunction.
    catch
        % In case of trouble make the laplacian definite
        [Steklov_evecs, Steklov_evals] = eigs(WW - 1e-8*speye(size(WW)), SS, Steklov_settings.Steklov_num_eigs + 1, 'sm');
        warning('The Steklov boundary is likely bad.')
    end
    
%     if min(min(evals)) < 0 %Correct a slightly negative lowest eigenvalue.
%        evals(1,1) = 0; 
%     end

    Steklov_evals(1,:) = []; %REMOVE THE CONSTANT EIGENFUNCTION AND EIGENVALUE.
    Steklov_evals(:,1) = [];
    
    Steklov_evecs(:,1) = [];
    
    prods = Steklov_evecs'*SS*Steklov_evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %L2(boundary) normalization, just in case (sometimes necessary) (WHY?)
    Steklov_evecs = Steklov_evecs * inv_pr;


%% Compute the Dirichlet-Steklov basis

DS_evecs = zeros(nv-length(landmarks), Steklov_settings.DS_num_eigs, length(landmarks) );
DS_evals = zeros(Steklov_settings.DS_num_eigs, Steklov_settings.DS_num_eigs, length(landmarks) );


for i = 1:length(landmarks) %Cycle over the landmarks. Each landmark gets the Steklov condition in turn.
            
    Dirichlet_boundary = vertcat( boundaries_segment{[1:(i-1) (i+1):length(landmarks)]} ); %All boundary indices, but that of the ith boundary, which is Steklov.
    
    WWW = WW; %Create the relevant operators.
    SSS = SS;    
    
    WWW(Dirichlet_boundary,:) = [];
    WWW(:,Dirichlet_boundary) = [];
    
    SSS(Dirichlet_boundary,:) = [];
    SSS(:,Dirichlet_boundary) = [];
        
    
    normSSS = sum(sum(SSS));
    SSS = SSS/normSSS; %Normalizes Steklov boundary length to 1. SUPER IMPORTANT
    
%     MassBoundary{i} = diag(diag( sum(S(boundaries_segment{i}, boundaries_segment{i}) , 1)) )/ normSS; %Lumped mass on boundary component.
    MassBoundary{i} = diag( sum(SS(boundaries_segment{i}, boundaries_segment{i}) , 1)) / normSSS; %Lumped mass on boundary component.
    
    try
        [pre_evecs, pre_evals] = eigs(WWW, SSS, Steklov_settings.DS_num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        [pre_evecs, pre_evals] = eigs(WWW - 1e-8*speye(size(WWW)), SSS, Steklov_settings.DS_num_eigs, 'sm');
        warning('The Steklov boundary is likely bad.')
    end

    prods = pre_evecs'*SSS*pre_evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %Normalization, just in case (sometimes necessary) (WHY?)
    
    pre_evecs = pre_evecs * inv_pr;
    
    sign_of_first_ef = sign( sum( sign( pre_evecs(:,1) ) ) ); % Fix a positive sign for the first eigenfunction.
    pre_evecs(:,1) = sign_of_first_ef*pre_evecs(:,1);
    
    DS_evecs(:,:,i) = BoundaryReinsertBatch(pre_evecs, Dirichlet_boundary, zeros(size(Dirichlet_boundary)));
    
    DS_evals(:,:,i) = pre_evals; %Normalization by the largest value. TODO: test other normalizations.
%     DS_evals(:,:,i) = pre_evals/max(max(pre_evals)); %Normalization by the largest value. TODO: test other normalizations.
  
    
    
end
    
    
    
%% Determine the interior vertices -- Initial interior vertices

Interior = 1:Shape.SHAPE.nv;
Interior(landmarks) = [];


%% Determine the refined interior vertices (not an output)

segment = 1:Shape_refined.SHAPE.nv;
segment(landmarks) = [];

%% Bases restricted to the boundary components.

for i = 1:length(landmarks) % Both bases restricted to the boundary.
    
    Shape_refined.STEKLOV.Steklov_evecs_bound{i} = Steklov_evecs(boundaries_segment{i}, : );  
    Shape_refined.STEKLOV.DS_evecs_bound{i} = DS_evecs(boundaries_segment{i}, : ,i);
%     size( Shape_refined.STEKLOV.evecs_bound{i})

end


%% Projection between the Steklov and DS bases (In L2(boundary) sense and normalization) (this is important)

P_DS_to_S = cell(length(landmarks),1); %Projection from the DS basis to the Steklov one.
                                       %The reverse projection is simply
                                       %the transpose of this one

for i = 1:length(landmarks)
    
    P_DS_to_S{i} = Shape_refined.STEKLOV.Steklov_evecs_bound{i}' * MassBoundary{i} * Shape_refined.STEKLOV.DS_evecs_bound{i};    
    
end







%% Combining the outputs into one neat structure
Shape_refined.STEKLOV.segment_TRIV = segment_TRIV;
% Shape_refined.STEKLOV.evecs = evecs;
Shape_refined.STEKLOV.Steklov_evecs = extendSegmentFunctionsIntegrals(Steklov_evecs, Shape_refined.STEKLOV.Steklov_evecs_bound, MassBoundary, segment, landmarks, Shape_refined.SHAPE.nv);
Shape_refined.STEKLOV.Steklov_evals = Steklov_evals/max(max(Steklov_evals));

Shape_refined.STEKLOV.DS_evecs = extendSegmentFunctionsIntegrals(DS_evecs, Shape_refined.STEKLOV.DS_evecs_bound, MassBoundary, segment, landmarks, Shape_refined.SHAPE.nv);
Shape_refined.STEKLOV.DS_evals = DS_evals;

Shape_refined.STEKLOV.boundary_edges_refined = boundary_edges_refined;
Shape_refined.STEKLOV.boundaries_segment = boundaries_segment;
Shape_refined.STEKLOV.MassBoundary = MassBoundary;
Shape_refined.STEKLOV.Interior = Interior;
Shape_refined.STEKLOV.boundaries = boundaries;
% Shape_refined.STEKLOV.normalTriangles = normalTriangles;

Shape_refined.STEKLOV.boundary_list = unique(boundary_edges_refined);
% Shape_refined.STEKLOV.harmonic_basis = ComputeRefinedCentralHarmonicBasis(Shape_refined, Shape_refined.STEKLOV.boundaries, landmarks); %Should no longer be necessary.

Shape_refined.STEKLOV.P_DS_to_S = P_DS_to_S;


%% Compute the Dirichlet LB basis and output it

[LB_evecs, LB_evals] = computeDirichletLB(Shape_refined.SHAPE, [landmarks(:); Shape_refined.STEKLOV.boundary_list(:)], 200);
LB_evecs = LB_evecs * diag( LB_evals.^-0.5 );

Shape_refined.STEKLOV.LB_evecs = LB_evecs;
Shape_refined.STEKLOV.LB_evals = LB_evals;



    
    
    


end