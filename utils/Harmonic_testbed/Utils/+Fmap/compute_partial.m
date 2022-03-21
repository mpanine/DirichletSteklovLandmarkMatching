function [ Fmap_info ] = compute_partial( Src, Tar, Fmap_settings, pair_settings, cache_settings, options )
    %COMPUTE_PARTIAL Computation of partial Fmaps

    numBasisFun_Src = pair_settings.basis_settings_Src.numBasisFun;
    numBasisFun_Tar = pair_settings.basis_settings_Tar.numBasisFun;
        
    if numBasisFun_Src ~= numBasisFun_Tar
        throw(MException('FMAP_COMPUTEPARTIALERROR:numBasisFunDifferent',...
                        'The bases on Source and Target must have the same size to use the partial Fmap computation.'));
    end
    options = struct;
    options.n_eigen = numBasisFun_Src;       % no. of eigenfunctions to use
    options.max_iters = pair_settings.max_iter;      % max no. iterations for the alternating optimization
    options.icp_max_iters = pair_settings.icp_max_iter; % max no. iterations for intrinsic ICP
    options.refine_fps = pair_settings.refine_fps;    % no. of fps for the final refinement
    options.fps_variance = pair_settings.fps_variance; % variance of gaussians for the final refinement

    % WARNING: the parameters below were chosen for the TOSCA dataset, and are
    % scale-dependent. To use them as-is, make sure your shapes have surface
    % area around 1.5 - 2.0 * 10^4.
    options.tv_sigma = Fmap_settings.energy_settings.tv_sigma;     % variance of indicator func. in TV term
    options.tv_mean = Fmap_settings.energy_settings.tv_mean;      % mean =
    

    options.mu1 = Fmap_settings.energy_settings.weight_slantedDiagonalMask;   % slanted-diagonal mask
    options.mu2 = Fmap_settings.energy_settings.weight_subOrthogonality; % sub-orthogonality
    options.mu3 = Fmap_settings.energy_settings.weight_area;   % area
    options.mu4 = Fmap_settings.energy_settings.weight_regularity; % regularity
    
    
    [C, ~, matches] = match_part_to_whole(Src, Tar, Src.DESCRIPTORS.fct,...
                                        Tar.DESCRIPTORS.fct,pair_settings.C_init, options);
    
    [C, ~, matches] = Fmap.refine_partial(Src, Tar, C, matches, options);
    
    Fmap_info.C = C;
    Fmap_info.matches = matches;
    Fmap_info.settings = Fmap_settings;
    Fmap_info.pair_settings = pair_settings;
    
end

