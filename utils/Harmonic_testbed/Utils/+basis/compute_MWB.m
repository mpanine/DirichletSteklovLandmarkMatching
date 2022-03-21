function [BASIS] = compute_MWB(S, MWB_settings)
    %COMPUTE_MWB Computes the manifold wavelet basis
    
    if ~isfield(S, 'sampleIndices')
        throw(MException('COMPUTE_LDBERROR:potentialNotDefined',...
                         sprintf('sampleIndices needs to be defined on the shape %s to allow the computation of the MWB.',...
                                    S.name)));
    end
    
    numWavelets = MWB_settings.parameters.numWaveletsPerScale;
    numLevels = length(MWB_settings.parameters.timeScales);
    type = MWB_settings.parameters.type;
    
    % Rescaling the timescales, so that they are independant of the shape
    % size
%     time_factor = mean(mean(S.A))/mean(mean(abs(S.W)))
%     max_geo_dist = max(max(S.Gamma(~isinf(S.Gamma))));
%     S.Gamma(isinf(S.Gamma)) = max_geo_dist;
%     mean_node_neigh_dist = cellfun(@(i_neighs,i)mean(S.Gamma(i,i_neighs)),S.vtx_neigh,mat2cell([1:S.nv]',ones(S.nv,1),1));
%     mean_node_neigh_dist(isnan(mean_node_neigh_dist)) = 0;
%     time_factor = mean(mean_node_neigh_dist);

    min_node_neigh_dist = cellfun(@(i_neighs,i) min(S.Gamma(i,i_neighs)),S.vtx_neigh,num2cell(1:S.nv,1)');
    min_node_neigh_dist(isnan(min_node_neigh_dist)) = 0;

    min_edge_length = min(min_node_neigh_dist);
    
    ts = MWB_settings.parameters.timeScales*min_edge_length;
    
    numBasisFun = MWB_settings.numBasisFun;
    
    samples = S.sampleIndices(1:numWavelets);
    
    orthonormalization = MWB_settings.parameters.orthonormalization;
    
    
    scale_fun = {};
    scale_fun_h = {};
    wavelet_fun = {};

    h_0 = 1e-3;
    h=h_0*ts(1);

    diracs = speye(S.nv);
    diracs = diracs(:,S.sampleIndices(1:numWavelets));

%     previous_phi = (eye(S.nv)+(ts(1)/2)*S.W)\(S.A*diracs);
    previous_phi = (S.A+(ts(1)/2)*S.W)\(S.A*diracs);
    previous_phi = bsxfun(@rdivide,previous_phi,sum(S.A*abs(previous_phi)));
%     previous_phi = [diracs,ones(S.nv,1)];

    for i=1:numLevels
        h=h_0*ts(i);
        % applying T^2 to phi_{j-1}
        phi_tmp = (S.A+ts(i)*S.W)\(S.A*previous_phi);
        phi_tmp_h = (S.A+(ts(i)+h)*S.W)\(S.A*previous_phi);
%         phi_tmp_2h = (S.A+(ts(i)+2*h)*S.W)\(S.A*previous_phi);
%         phi_tmp_3h = (S.A+(ts(i)+3*h)*S.W)\(S.A*previous_phi);

        phi_tmp = bsxfun(@rdivide,phi_tmp,sum(S.A*abs(phi_tmp)));
        phi_tmp_h = bsxfun(@rdivide,phi_tmp_h,sum(S.A*abs(phi_tmp_h)));
%         phi_tmp_2h = bsxfun(@rdivide,phi_tmp_2h,sum(S.A*abs(phi_tmp_2h)));
%         phi_tmp_3h = bsxfun(@rdivide,phi_tmp_3h,sum(S.A*abs(phi_tmp_3h)));

        scale_fun{end+1} = phi_tmp;
        scale_fun_h{end+1} = phi_tmp_h;

%         wavelet_fun{end+1} = (-phi_tmp_3h+5*phi_tmp_2h-7*phi_tmp_h+3*phi_tmp)/(2*h^2); % old formulation

%         wavelet_fun{end+1} = (2*phi_tmp-5*phi_tmp_h+4*phi_tmp_2h-phi_tmp_3h)/h^2; % correct version of the old formulation

%         wavelet_fun{end+1} = (phi_tmp-2*phi_tmp_h+1*phi_tmp_2h)/h^2;
        wavelet_fun{end+1} = phi_tmp-phi_tmp_h;

        previous_phi = scale_fun{end};
    end


    basis = [];
    scales = [];

    switch type
        case 'mexican hat'
            % MEXICAN HAT wavelets
            for i=1:numLevels
                basis = [basis,wavelet_fun{i}];
                scales = [scales,scale_fun{i}];
            end
        case 'haar'
            % HAAR wavelets
            scales2 = [];
            mean_val = mean(scale_fun{1},1);
            for i=1:numLevels
                fun1 = scale_fun{i};
                fun1(fun1>=mean_val(1)) = 1;
                fun1(fun1<mean_val(1)) = 0;
                fun2 = scale_fun_h{i};
                fun2(fun2>=mean_val(1)) = 1;
                fun2(fun2<mean_val(1)) = 0;
                scales = [scales,fun1];
                scales2 = [scales2,fun2];
            end
            basis = 2*scales-scales2;
    end
    
    if orthonormalization
        basis = weighted_orthonormal_basis_noPrint(basis, S.A);
    end

    % Removing the NaN values
    basis(isnan(basis)) = 1;
    % Computing singular values for the basis as "eigenvalues"
    try
        [~,singular_mat,~] = svd(basis);
    catch error
        warning('NaN or Inf values during the svd computation of the MWB.');
        singular_mat = NaN*ones(size(basis,2));
    end
    
    basis = full(basis);
    
    BASIS = struct;
    BASIS.basis = basis;
    BASIS.basis_inverse = basis'*S.A;
    BASIS.scales = scales;
    BASIS.eigenvalues = diag(singular_mat);
    
end

