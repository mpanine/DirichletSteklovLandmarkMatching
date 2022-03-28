function [Shape_refined, segment_TRIV, Steklov_evecs, Steklov_evals, boundary_edges_old, boundary_old] = ComputePrincipledSteklovAll(Shape, landmarks, Steklov_settings)
%Computes the both the Steklov and Dirichlet-Steklov bases for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark



Shape_refined = Shape;

[Shape_refined.SHAPE, boundaries] = refineCircleAroundCenter(Shape.SHAPE, landmarks, Steklov_settings.landmarks_radii, Steklov_settings.numsubdivisions); %Insert circles at landmarks.

nv = size(Shape_refined.SHAPE.surface.VERT,1); %Just to make sure the size is correctly set.
Shape_refined.SHAPE.nv = nv;
Shape_refined.SHAPE.nv_unrefined = Shape.SHAPE.nv; % Original number of vertices.


boundary_old = vertcat(boundaries{:});%cell2mat(boundary); %List form of the boundary indices: loses information about separate segments. TODO: output all.

[WW, SS, segment_TRIV, boundaries_segment, boundary_edges_refined] = splitMeshLandmarkSteklov_clean(Shape_refined, landmarks, boundaries);


%% Normalize mass matrices and compute MassBoundary 

MassBoundary = cell(1,length(landmarks));

for i = 1:length(landmarks) %Cycle over the landmarks. Each landmark gets the Steklov condition in turn.            
  
    SSS = diag( sum(SS(boundaries_segment{i}, boundaries_segment{i}) , 1));
    normSSS = sum(sum(SSS));

    MassBoundary{i} = SSS;%/normSSS;
    
    SS(boundaries_segment{i}, boundaries_segment{i}) = SS(boundaries_segment{i}, boundaries_segment{i});% / normSSS; %Normalizes each boundary component.
    
end

Shape_refined.STEKLOV.S = operatorBoundaryReinsert(SS, landmarks); %Outputs the normalized Steklov mass matrix.

% total_boundary_length = sum(sum(SS))

%% Determine the refined interior vertices

segment = 1:Shape_refined.SHAPE.nv;
segment(landmarks) = [];

%% Determine the interior vertices -- Initial interior vertices

Interior = 1:Shape.SHAPE.nv;
Interior(landmarks) = [];

%% Compute the Dirichlet-Steklov basis



DS_evecs = zeros(nv-length(landmarks), Steklov_settings.DS_num_eigs, length(landmarks) );
DS_evals = zeros(Steklov_settings.DS_num_eigs, Steklov_settings.DS_num_eigs, length(landmarks) );
DS_sqrt_inv_evals = DS_evals; % inverse square root of evals.


for i = 1:length(landmarks) %Cycle over the landmarks. Each landmark gets the Steklov condition in turn.
            
    Dirichlet_boundary = vertcat( boundaries_segment{[1:(i-1) (i+1):length(landmarks)]} ); %All boundary indices, but that of the ith boundary, which is Steklov.
    
    WWW = WW; %Create the relevant operators.
    SSS = SS;    
    
    WWW(Dirichlet_boundary,:) = [];
    WWW(:,Dirichlet_boundary) = [];
    
    SSS(Dirichlet_boundary,:) = [];
    SSS(:,Dirichlet_boundary) = [];
        
    
    normSSS = sum(sum(SSS));
    SSS = SSS;%/normSSS; %Normalizes Steklov boundary length to 1. SUPER IMPORTANT
    
    MassBoundary{i} = diag( sum(SS(boundaries_segment{i}, boundaries_segment{i}) , 1));% / normSSS; %Lumped mass on boundary component.
    
    try
        [pre_evecs, pre_evals] = eigs(WWW, SSS, Steklov_settings.DS_num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        try
            [pre_evecs, pre_evals] = eigs(WWW - 1e-8*speye(size(WWW)), SSS, Steklov_settings.DS_num_eigs, 'sm');
        catch ME
            fprintf('Error at boundary %d of shape %s\n',i,Shape.SHAPE.name);
            rethrow(ME);
        end
        warning('The Steklov boundary is likely bad.')
    end
    pre_evals = diag(pre_evals);
    [pre_evals, pre_order] = sort(abs(pre_evals),'ascend');
    pre_evecs = pre_evecs(:,pre_order);
    pre_evals = diag(pre_evals);
    
    prods = pre_evecs'*SSS*pre_evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %L2 Normalization, just in case (sometimes necessary)
    
    pre_evecs = pre_evecs * inv_pr;
    
    sign_of_first_ef = sign( sum( sign( pre_evecs(:,1) ) ) ); % Fix a positive sign for the first eigenfunction.
    pre_evecs(:,1) = sign_of_first_ef*pre_evecs(:,1);
    
    DS_evecs(:,:,i) = BoundaryReinsertBatch(pre_evecs, Dirichlet_boundary, zeros(size(Dirichlet_boundary))) * diag( diag(pre_evals).^-.5 ); %Reinsert boundaries and W-normalize.
    
    DS_evals(:,:,i) = pre_evals; 
    DS_sqrt_inv_evals(:,:,i) = diag( diag(pre_evals).^-.5 );
    
end
    

    

%% Compute the Landmark harmonics and their normal derivatives



landmark_harmonics = ComputeRefinedCentralHarmonicBasis(Shape_refined, boundaries, landmarks); % This is defined on the full shape
landmark_harmonics_ND = cell(length(landmarks), 1);


pre_lh_ND = WW * landmark_harmonics(segment,:);


for i = 1:length(landmarks)
    
%     unproj_landmark_harmonics_ND{i} = pre_lh_ND(boundaries_segment{i}, :);
    landmark_harmonics_ND{i} =  DS_evecs(boundaries_segment{i},:,i) * sqrt(DS_evals(:,:,i)) * DS_evecs(boundaries_segment{i},:,i)' * MassBoundary{i} * pre_lh_ND(boundaries_segment{i}, :);
    
end





%% Bases restricted to the boundary components.

for i = 1:length(landmarks) % Both bases restricted to the boundary.
    
%     Shape_refined.STEKLOV.Steklov_evecs_bound{i} = Steklov_evecs(boundaries_segment{i}, : );  
    Shape_refined.STEKLOV.DS_evecs_bound{i} = DS_evecs(boundaries_segment{i}, : ,i);
%     size( Shape_refined.STEKLOV.evecs_bound{i})

end




%% Compute the Dirichlet LB basis

Shape_refined.STEKLOV.boundary_list = unique(boundary_edges_refined); %Important output and useful here.

[LB_evecs, LB_evals] = computeDirichletLB(Shape_refined.SHAPE, [landmarks(:); Shape_refined.STEKLOV.boundary_list(:)], Steklov_settings.num_LB_eigs,Steklov_settings.num_LB_margin,Steklov_settings.eig_threshold);
LB_evecs = LB_evecs * diag( LB_evals.^-0.5 );

Shape_refined.STEKLOV.laplaceBasis = struct;
Shape_refined.STEKLOV.laplaceBasis.basis = LB_evecs;
Shape_refined.STEKLOV.laplaceBasis.basis_inverse = LB_evecs';
Shape_refined.STEKLOV.laplaceBasis.eigenvalues = LB_evals;

LB_evecs(landmarks,:) = []; %REMOVE THE LANDMARKS




%% Combine and output the full basis (first DS in landmark order, then LB)


FullBasis = zeros(nv - length(landmarks), length(landmarks)*Steklov_settings.DS_num_eigs + Steklov_settings.num_LB_eigs);

for i = 1:length(landmarks)
    
    ind_start = (i-1)*Steklov_settings.DS_num_eigs + 1;
    ind_end = i*Steklov_settings.DS_num_eigs;
    
    FullBasis(:, ind_start:ind_end) = DS_evecs(:,:,i);
    
    
end

FullBasis(:, (length(landmarks)*Steklov_settings.DS_num_eigs + 1):end) = LB_evecs;

Shape_refined.STEKLOV.FullBasisProduct = FullBasis' * WW * FullBasis; % Matrix of W-inner products.

Shape_refined.STEKLOV.FullBasis_segment = FullBasis;
Shape_refined.STEKLOV.FullBasis = extendSegmentFunctionsIntegrals_principled(FullBasis, boundaries_segment, MassBoundary, segment, landmarks, Shape_refined.SHAPE.nv);






%% Combining the outputs into one neat structure
Shape_refined.STEKLOV.segment_TRIV = segment_TRIV;

Shape_refined.STEKLOV.DS_evecs = extendSegmentFunctionsIntegrals(DS_evecs, Shape_refined.STEKLOV.DS_evecs_bound, MassBoundary, segment, landmarks, Shape_refined.SHAPE.nv);
Shape_refined.STEKLOV.DS_evals = DS_evals;
Shape_refined.STEKLOV.DS_sqrt_inv_evals = DS_sqrt_inv_evals;

Shape_refined.STEKLOV.boundary_edges_refined = boundary_edges_refined;
Shape_refined.STEKLOV.boundaries_segment = boundaries_segment;
Shape_refined.STEKLOV.MassBoundary = MassBoundary;
Shape_refined.STEKLOV.Interior = Interior;
Shape_refined.STEKLOV.boundaries = boundaries;


Shape_refined.STEKLOV.S_segment = SS;
Shape_refined.STEKLOV.W_segment = WW;
Shape_refined.STEKLOV.segment = segment;


Shape_refined.STEKLOV.landmark_harmonics = landmark_harmonics;
Shape_refined.STEKLOV.landmark_harmonics_ND = landmark_harmonics_ND;

% Shape_refined.STEKLOV.LB_evecs = LB_evecs;
% Shape_refined.STEKLOV.LB_evals = LB_evals;

Shape_refined.STEKLOV.laplaceBasis_segment = struct;
Shape_refined.STEKLOV.laplaceBasis_segment.basis = LB_evecs;
Shape_refined.STEKLOV.laplaceBasis_segment.basis_inverse = LB_evecs';
Shape_refined.STEKLOV.laplaceBasis_segment.eigenvalues = LB_evals;
end
