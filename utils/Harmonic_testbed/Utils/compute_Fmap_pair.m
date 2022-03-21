function [ C_Fmap_PreRef, ...
           C_Fmap_ZO, ...
           pF_Fmap_PreRef, ...
           pF_Fmap_ZO, time_INI, time_ZO] = compute_Fmap_pair( Src, Tar, ...
                                    dataset,...
                                    descriptor_settings, ...
                                    Fmap_settings,...
                                    p2p_settings,...
                                    Fmap_ref_params,...
                                    cache_settings, ...
                                    basis_settings_Fmap_ini,...
                                    basis_settings, ...
                                    energy_settings)
    %COMPUTE_FMAP_PAIR Computes the Fmap given a Source-Target pair.
    Src_Obj.SHAPE = Src;
    Tar_Obj.SHAPE = Tar;
    
    start_time = tic;
    Src_Obj.DESCRIPTORS = descriptor.compute(Src_Obj.SHAPE, descriptor_settings, cache_settings);
    Tar_Obj.DESCRIPTORS = descriptor.compute(Tar_Obj.SHAPE, descriptor_settings, cache_settings);

    Src_Obj.BASIS_ZoomOut = basis.compute(Src_Obj.SHAPE, basis_settings, cache_settings);
    Tar_Obj.BASIS_ZoomOut = basis.compute(Tar_Obj.SHAPE, basis_settings, cache_settings);
    
    Src_Obj.BASIS = struct;
    Tar_Obj.BASIS = struct;
    
    k = basis_settings_Fmap_ini.numBasisFun;
    
    Src_Obj.BASIS.basis = Src_Obj.BASIS_ZoomOut.basis(:,1:k);
    Src_Obj.BASIS.basis_inverse = Src_Obj.BASIS_ZoomOut.basis_inverse(1:k,:);
    Src_Obj.BASIS.eigenvalues = Src_Obj.BASIS_ZoomOut.eigenvalues(1:k);
    Src_Obj.BASIS.settings = basis_settings_Fmap_ini;
    Tar_Obj.BASIS.basis = Tar_Obj.BASIS_ZoomOut.basis(:,1:k);
    Tar_Obj.BASIS.basis_inverse = Tar_Obj.BASIS_ZoomOut.basis_inverse(1:k,:);
    Tar_Obj.BASIS.eigenvalues = Tar_Obj.BASIS_ZoomOut.eigenvalues(1:k);
    Tar_Obj.BASIS.settings = basis_settings_Fmap_ini;
    
    Src_Obj.basis_settings = basis_settings_Fmap_ini;
    Tar_Obj.basis_settings = basis_settings_Fmap_ini;
    
    [Energies, dEnergies] = energy.compute(Src_Obj, Tar_Obj, energy_settings);
    [objectiveFunction, objectiveGradient] = ...
        energy.compute_objectiveFunction(Energies, dEnergies);

    pair_settings = shape.PairSettings(Src_Obj, Tar_Obj, dataset,...
                                       objectiveFunction, objectiveGradient);

    Fmap_info = Fmap.compute(Src_Obj, Tar_Obj, Fmap_settings, pair_settings, cache_settings);
    C_Fmap_PreRef = Fmap_info.C;
    pF_Fmap_PreRef = p2pMap.compute(Src_Obj.BASIS,Tar_Obj.BASIS,C_Fmap_PreRef,p2p_settings);
    time_INI = toc(start_time);
    
%     pF_Fmap_PreRef = Fmap_info.matches;
    start_time = tic;
    [C_Fmap_ZO, pF_Fmap_ZO] = Fmap.refine(Src_Obj, Tar_Obj, C_Fmap_PreRef, 'type', 'ZoomOut', 'parameters', Fmap_ref_params);
    time_ZO = toc(start_time)+time_INI;
    
end
