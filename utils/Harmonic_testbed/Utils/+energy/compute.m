function [ varargout ] = compute( Src, Tar, energy_settings )
    %COMPUTE Computes the energy terms, depending on the type of
    %energy_settings provided.
    
    fprintf('Computing Energy with settings %s...',class(energy_settings)); tic;
    
    switch class(energy_settings)
        case 'energy.Full_Settings'
            [Energies, dEnergies] = energy.compute_full(Src, Tar, energy_settings);
            energy_terms = {Energies, dEnergies};
        case 'energy.Partial_Settings'
            [ Energy_F, dEnergy_F, Energy_v, dEnergy_v, x0_v, F_init, est_rank ] = energy.compute_partial(Src, Tar, energy_settings);
            energy_terms = {Energy_F, dEnergy_F, Energy_v, dEnergy_v, x0_v, F_init, est_rank};
        otherwise
            throw(MException('ENERGYCOMPUTEERROR:energy_settingsUnknown',...
                             'The energy_settings object provided is unknown.'));
    end;
    varargout = cell(size(energy_terms));
    for i=1:size(energy_terms, 2)
        varargout{i} = energy_terms{i};
    end
    
    t =toc; fprintf('done:%.4fs\n',t);
end

