function [ time_INI,time_ZO ] = compute_Fmap_pair_2_files( Src, Tar, ...
                                    dataset,...
                                    i,...
                                    descriptor_settings, ...
                                    Fmap_settings,...
                                    p2p_settings,...
                                    Fmap_ref_params,...
                                    cache_settings, ...
                                    basis_settings_Fmap_INI,...
                                    basis_settings, ...
                                    energy_settings, ...
                                    saveDir, fun2useStr )
    fun2use = str2func(fun2useStr);
    %COMPUTE_FMAP_PAIR Computes the Fmap given a Source-Target pair. The
    %output can be saved to a file.
    
    % /!\ Fmap computed from Tar to Src so that the pMap goes from Src to
    % Target.
    cache_settings.enabled = false;
    [ C_Fmap_INI, ...
       C_Fmap_ZO, ...
       pF_Fmap_INI, ...
       pF_Fmap_ZO,time_INI,time_ZO] = fun2use( Tar, Src, ...
                                dataset,...
                                descriptor_settings, ...
                                Fmap_settings,...
                                p2p_settings,...
                                Fmap_ref_params,...
                                cache_settings, ...
                                basis_settings_Fmap_INI, ...
                                basis_settings, ...
                                energy_settings);
    cache_settings.enabled = true;
    
    save(sprintf('%s/%d_%s-%s_pF_Fmap_INI.mat',char(saveDir), i, Src.name, Tar.name),'pF_Fmap_INI');
    save(sprintf('%s/%d_%s-%s_pF_Fmap_ZO.mat',char(saveDir), i, Src.name, Tar.name),'pF_Fmap_ZO');
    save(sprintf('%s/%d_%s-%s_C_Fmap_INI.mat',char(saveDir), i, Src.name, Tar.name),'C_Fmap_INI');
    save(sprintf('%s/%d_%s-%s_C_Fmap_ZO.mat',char(saveDir), i, Src.name, Tar.name),'C_Fmap_ZO');
    
end

