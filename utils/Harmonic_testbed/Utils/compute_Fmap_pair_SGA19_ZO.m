function [ C_Fmap_PreRef, ...
           C_Fmap_ZO, ...
           pF_Fmap_PreRef, ...
           pF_Fmap_ZO, time_INI, time_ZO] = compute_Fmap_pair_SGA19_ZO( Src, Tar, ...
                                    dataset,...
                                    descriptor_settings, ...
                                    Fmap_settings,...
                                    p2p_settings,...
                                    Fmap_ref_params,...
                                    cache_settings, ...
                                    basis_settings_Fmap_ini,...
                                    basis_settings, ...
                                    energy_settings)
    %COMPUTE_FMAP_PAIR_SGA19_ZO Computation of Fmap from Src to Tar using
    %the code from SGA19 zoomOut.
    numLmks = length(descriptor_settings.landmarks);
    numEigsDsc = descriptor_settings.parameters.numEigs;
    t = descriptor_settings.parameters.numTimes;
    numSkip = descriptor_settings.parameters.numSkip;
    a = energy_settings.weight_descriptorPreservation;
    b = energy_settings.weight_descriptorCommutativity;
    c = energy_settings.weight_LBCommutativity;
    k_init = basis_settings_Fmap_ini.numBasisFun;
    k = basis_settings.numBasisFun;
    
    start_time = tic;
    BASIS = basis.compute(Src, basis_settings, cache_settings);
    B1 = struct;
    B1.evecs = BASIS.basis;
    B1.evals = BASIS.eigenvalues;
    Src.evecs = B1.evecs;
    Src.evals = B1.evals;
    BASIS = basis.compute(Tar, basis_settings, cache_settings);
    B2 = struct;
    B2.evecs = BASIS.basis;
    B2.evals = BASIS.eigenvalues;
    Tar.evecs = B2.evecs;
    Tar.evals = B2.evals;
    
    fct_Src = fMAP.compute_descriptors_with_landmarks(Src,numEigsDsc,Src.landmarkIndices(1:numLmks),t,numSkip);
    fct_Tar = fMAP.compute_descriptors_with_landmarks(Tar,numEigsDsc,Tar.landmarkIndices(1:numLmks),t,numSkip);
    para = struct;
    para.type_orient = 'direct';
    para.a = a;para.b = b;para.c = c;  
    
    para.fMap_size = [k_init, k_init];
    C_Fmap_PreRef = fMAP.compute_fMap_from_descriptors(Src, Tar, fct_Src, fct_Tar, para);
    pF_Fmap_PreRef = fMAP.fMap2pMap(B1,B2,C_Fmap_PreRef);

    time_INI = toc(start_time);
    
    % apply zoomOut
    para.k_init = k_init;
    para.k_step = 1;
    para.k_final = k;
    start_time = tic;
    [pF_Fmap_ZO, C_Fmap_ZO, ~, ~] = zoomOut_refine(Tar.evecs, Src.evecs, pF_Fmap_PreRef, para);
    time_ZO = toc(start_time)+time_INI;
    fprintf('ZoomOut runtime: %.2f sec\n',t)
end

