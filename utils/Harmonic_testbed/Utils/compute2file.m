function [ output_name ] = compute2file( singleBasisEvaluation_settings, directory, name )
    %COMPUTE2FILE Same as compute, but writes content on disk as it goes,
    %based either on a name specified or a default name, generated based on
    %the date of execution.
    
    if nargin < 3
        name = [singleBasisEvaluation_settings.dataset, ' ', strrep(datestr(datetime), ':','-')];
    end
    
    output_name = name;
    
    saveDir = [directory,name]
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    
    assert(isa(singleBasisEvaluation_settings, 'singleBasisEvaluation.Settings'),...
        sprintf('singleBasisEvaluation_settings has to be a singleBasisEvaluation.Settings object got %s object instead.',...
                                                    class(singleBasisEvaluation_settings)));

    save([saveDir,'/settings.mat'],'singleBasisEvaluation_settings');
    
    
    
    nDescr = length(singleBasisEvaluation_settings.descriptorsIdx);
    nEval = length(singleBasisEvaluation_settings.evaluationIdx);

    for i=1:length(singleBasisEvaluation_settings.shapePairList)
        fprintf(sprintf('Iteration %d/%d\n',i,length(singleBasisEvaluation_settings.shapePairList)));
        pair_str = singleBasisEvaluation_settings.shapePairList(i);
        [names, ~] = split(pair_str);
        shape_name_Src = char(names(1));
        shape_name_Tar = char(names(2));

        switch singleBasisEvaluation_settings.dataset
            case {'Faust_remeshed', 'SCAPE_remeshed', 'TOSCA_Isometric', 'TOSCA_nonIsometric'}
                mesh_dir = sprintf('Data/%s/vtx_5k/',singleBasisEvaluation_settings.dataset);
            case 'Shrec19_remeshed'
                mesh_dir = sprintf('Data/%s/vtx_5k/remeshed_rescaled/',singleBasisEvaluation_settings.dataset);
            otherwise
                mesh_dir = sprintf('Data/%s/',singleBasisEvaluation_settings.dataset);
        end

        % Settings
        shape_settings = singleBasisEvaluation_settings.shape_settings;
        cache_settings = singleBasisEvaluation_settings.cache_settings;
        useCache = cache_settings.enabled;

        descriptor_settings_Src = singleBasisEvaluation_settings.descriptor_settings;
        descriptor_settings_Tar = singleBasisEvaluation_settings.descriptor_settings;

        switch singleBasisEvaluation_settings.dataset
            case 'Shrec16_Partial/holes'
                split_ls = split(mesh_dir, '/');
                dataDir =  char(split_ls(1));
                lmkDir =  char(split_ls(2));
                landmarks_Src = load([dataDir,'/',lmkDir,'/Landmarks/','holes/',shape_name_Src,'.mat']);
                landmarks_Tar = load([dataDir,'/',lmkDir,'/Landmarks/','holes/',shape_name_Tar,'.mat']);
            case 'Shrec16_Partial/cuts'
                split_ls = split(mesh_dir, '/');
                dataDir =  char(split_ls(1));
                lmkDir = char(split_ls(2));
                landmarks_Src = load([dataDir,'/',lmkDir,'/Landmarks/','cuts/',shape_name_Src,'.mat']);
                landmarks_Tar = load([dataDir,'/',lmkDir,'/Landmarks/','cuts/',shape_name_Tar,'.mat']);
            case {'Faust_remeshed', 'SCAPE_remeshed', 'TOSCA_Isometric', 'TOSCA_nonIsometric'}
                dataDir = sprintf('Data/%s/',singleBasisEvaluation_settings.dataset);
                landmarks_Src = load([dataDir,'Landmarks/',shape_name_Src,'.mat']);
                landmarks_Tar = load([dataDir,'Landmarks/',shape_name_Tar,'.mat']);
            case {'Shrec19','Shrec19_remeshed'}
                dataDir = sprintf('Data/%s/',singleBasisEvaluation_settings.dataset);
                landmarks_Src = load([dataDir,'Landmarks/',shape_name_Src,'.mat']);
                landmarks_Tar = load([dataDir,'Landmarks/',shape_name_Src,'_',shape_name_Tar,'.mat']);
            otherwise
                landmarks_Src = load([mesh_dir,'Landmarks/',shape_name_Src,'.mat']);
                landmarks_Tar = load([mesh_dir,'Landmarks/',shape_name_Tar,'.mat']);
        end

        if strcmp(singleBasisEvaluation_settings.metricType,'Distortion')
            if length(landmarks_Src.landmarks) < nDescr
                throw(MException('COMPUTESINGLEBASISEVALUATION:notEnoughLandmarks_Src',...
                    'Not enough landmarks defined on Source shape to allow computation: the number of landmarks defined must be greater than singleBasisEvaluation.Settings.numLandmarksDescriptors for the chosen metric.'));
            end
            if length(landmarks_Tar.landmarks) < nDescr
                throw(MException('COMPUTESINGLEBASISEVALUATION:notEnoughLandmarks_Tar',...
                    'Not enough landmarks defined on Tar shape to allow computation: the number of landmarks defined must be greater than singleBasisEvaluation.Settings.numLandmarksDescriptors for the chosen metric.'));
            end
        else % the other metrics require to have enough evaluation landmarks (at least 10)
            if length(landmarks_Src.landmarks) < (nDescr+10)
                throw(MException('COMPUTESINGLEBASISEVALUATION:notEnoughLandmarks_Src',...
                    'Not enough landmarks defined on Source shape to allow computation: the number of landmarks defined must be greater than singleBasisEvaluation.Settings.numLandmarksDescriptors + 10 for the chosen metric.'));
            end
            if length(landmarks_Tar.landmarks) < (nDescr+10)
                throw(MException('COMPUTESINGLEBASISEVALUATION:notEnoughLandmarks_Tar',...
                    'Not enough landmarks defined on Tar shape to allow computation: the number of landmarks defined must be greater than singleBasisEvaluation.Settings.numLandmarksDescriptors + 10 for the chosen metric.'));
            end
        end
        if nDescr > 0
            descriptor_settings_Src.landmarks = landmarks_Src.landmarks(singleBasisEvaluation_settings.descriptorsIdx);
            descriptor_settings_Tar.landmarks = landmarks_Tar.landmarks(singleBasisEvaluation_settings.descriptorsIdx);
            % VERY IMPORTANT: activate the "haveLmks" flag of both
            % descriptor.Settings objects so that landmark-dependant descriptor
            % functions are computed.
            descriptor_settings_Src.haveLmks = true;
            descriptor_settings_Tar.haveLmks = true;
        end;

        if nEval > 0
            if length(landmarks_Src.landmarks) < (nDescr+nEval)
                throw(MException('COMPUTESINGLEBASISEVALUATION:notEnoughLandmarks_Src',...
                    'Not enough landmarks defined on Source shape to allow computation: the number of landmarks defined must be greater than singleBasisEvaluation.Settings.numLandmarksDescriptors + singleBasisEvaluation.Settings.numLandmarksEvaluation for the chosen metric.'));
            end
            if length(landmarks_Tar.landmarks) < (nDescr+nEval)
                throw(MException('COMPUTESINGLEBASISEVALUATION:notEnoughLandmarks_Tar',...
                    'Not enough landmarks defined on Tar shape to allow computation: the number of landmarks defined must be greater than singleBasisEvaluation.Settings.numLandmarksDescriptors + singleBasisEvaluation.Settings.numLandmarksEvaluation for the chosen metric.'));
            end
            evalLmk_Src = landmarks_Src.landmarks(singleBasisEvaluation_settings.evaluationIdx);
            evalLmk_Tar = landmarks_Tar.landmarks(singleBasisEvaluation_settings.evaluationIdx);

        else
            % if nEval = 0, taking all remaining landmarks for
            % evaluation
            evalLmk_Src = landmarks_Src.landmarks;
            evalLmk_Src(singleBasisEvaluation_settings.descriptorsIdx) = [];
            evalLmk_Tar = landmarks_Tar.landmarks;
            evalLmk_Tar(singleBasisEvaluation_settings.descriptorsIdx) = [];
        end

        basis_settings = singleBasisEvaluation_settings.basis_settings;
        basis_settings_ICP = singleBasisEvaluation_settings.basis_settings_ICP;
        basis_settings_BCICP = singleBasisEvaluation_settings.basis_settings_BCICP;
        basis_settings_ZoomOut = singleBasisEvaluation_settings.basis_settings_ZoomOut;

        % ---
        % Defining Source and Target structures
        Src = struct;
        Tar = struct;
        Src.SHAPE = shape.compute([mesh_dir,shape_name_Src], ...
                            shape_settings, cache_settings);      
        Src.SHAPE.landmarkIndices = landmarks_Src.landmarks;

        Tar.SHAPE = shape.compute([mesh_dir,shape_name_Tar], ...
                            shape_settings, cache_settings);
        Tar.SHAPE.landmarkIndices = landmarks_Tar.landmarks;

        switch singleBasisEvaluation_settings.basisType
            case 'LMH'
                recon_error = sqrt((Src.SHAPE.surface.VERT(:,1)-Src.SHAPE.laplacianBasis.basis*(Src.SHAPE.laplacianBasis.basis'*Src.SHAPE.surface.VERT(:, 1))).^2+...
                (Src.SHAPE.surface.VERT(:,2)-Src.SHAPE.laplacianBasis.basis*(Src.SHAPE.laplacianBasis.basis'*Src.SHAPE.surface.VERT(:, 2))).^2+...
                (Src.SHAPE.surface.VERT(:,3)-Src.SHAPE.laplacianBasis.basis*(Src.SHAPE.laplacianBasis.basis'*Src.SHAPE.surface.VERT(:, 3))).^2);
                recon_error = recon_error./max(recon_error(:));

                Src.SHAPE.region = recon_error;

                clear recon_error;

                recon_error = sqrt((Tar.SHAPE.surface.VERT(:,1)-Tar.SHAPE.laplacianBasis.basis*(Tar.SHAPE.laplacianBasis.basis'*Tar.SHAPE.surface.VERT(:, 1))).^2+...
                (Tar.SHAPE.surface.VERT(:,2)-Tar.SHAPE.laplacianBasis.basis*(Tar.SHAPE.laplacianBasis.basis'*Tar.SHAPE.surface.VERT(:, 2))).^2+...
                (Tar.SHAPE.surface.VERT(:,3)-Tar.SHAPE.laplacianBasis.basis*(Tar.SHAPE.laplacianBasis.basis'*Tar.SHAPE.surface.VERT(:, 3))).^2);
                recon_error = recon_error./max(recon_error(:));

                Tar.SHAPE.region = recon_error;
            case 'HOB'
                LBO_sum = zeros(1, Src.SHAPE.nv);
                for idx_lb=2:5
                    LB = Src.SHAPE.laplacianBasis.basis(:, idx_lb)'.^2;
                    LBO_sum = LBO_sum + ((LB - min(LB))/(max(LB) - min(LB)));
                end
                Src.SHAPE.potential = LBO_sum;

                clear LBO_sum;

                LBO_sum = zeros(1, Tar.SHAPE.nv);
                for idx_lb=2:5
                    LB = Tar.SHAPE.laplacianBasis.basis(:, idx_lb)'.^2;
                    LBO_sum = LBO_sum + ((LB - min(LB))/(max(LB) - min(LB)));
                end
                Tar.SHAPE.potential = LBO_sum;
        end
        
        if isa(descriptor_settings_Src,'descriptor.Shot_Settings')
            descriptor_settings_Src.landmarks = 1:Src.SHAPE.nv;
        end
        if isa(descriptor_settings_Tar,'descriptor.Shot_Settings')
            descriptor_settings_Tar.landmarks = 1:Tar.SHAPE.nv;
        end

        if singleBasisEvaluation_settings.basisType == 'MWB'
            warning('Not using the descriptor settings, as the MWB is used. Using scaling functions instead.');
            numLevels = length(basis_settings.parameters.timeScales);
            log_ts = linspace(size(Src.SHAPE.landmarkIndices,2)/Src.SHAPE.nv, 1, numLevels);
            ts = exp(log_ts);
            basis_settings.parameters.timeScales = ts;
            Src.BASIS = basis.compute(Src.SHAPE, basis_settings, cache_settings);
            Tar.BASIS = basis.compute(Tar.SHAPE, basis_settings, cache_settings);
            Src.DESCRIPTORS.fct = Src.BASIS.scales;
            Tar.DESCRIPTORS.fct = Tar.BASIS.scales;
            
            size(Src.BASIS.basis)
            size(Src.BASIS.eigenvalues)
        else
            Src.BASIS = basis.compute(Src.SHAPE, basis_settings, cache_settings);
            Tar.BASIS = basis.compute(Tar.SHAPE, basis_settings, cache_settings);
            Src.DESCRIPTORS = descriptor.compute(Src.SHAPE, descriptor_settings_Src, cache_settings);
            Tar.DESCRIPTORS = descriptor.compute(Tar.SHAPE, descriptor_settings_Tar, cache_settings);
        end

        Src.basis_settings = basis_settings;
        Tar.basis_settings = basis_settings;
        % ---


        % Computing the Energies and their derivatives
        energy_settings = singleBasisEvaluation_settings.energy_settings;


        refinement_parameters = singleBasisEvaluation_settings.refinement_parameters;

        switch class(energy_settings)
            case 'energy.Full_Settings'
                [Energies, dEnergies] = energy.compute(Src, Tar, energy_settings);

                % Computing the associated objective function
                [objectiveFunction, objectiveGradient] = energy.compute_objectiveFunction(Energies, dEnergies);

                % ---
                % Fmap settings
                pair_settings = shape.PairSettings(Src, Tar, singleBasisEvaluation_settings.dataset,...
                                        objectiveFunction, objectiveGradient);

                Fmap_settings = Fmap.Settings(energy_settings,...
                                              'optimizer', singleBasisEvaluation_settings.optimizer);
                % ---

                % Computing the Fmap
                cache_settings.enabled = false;
                Fmap_info = Fmap.compute(Src, Tar, Fmap_settings, pair_settings, cache_settings);
                cache_settings.enabled = useCache;

                % p2pMap settings
                p2pMap_settings = p2pMap.Settings(Fmap_settings,...
                                                Src.basis_settings, Tar.basis_settings,...
                                                'type' ,'KNN');

                % raw p2pMap computation
                pF = p2pMap.compute(Src.BASIS, Tar.BASIS, Fmap_info.C, p2pMap_settings );

            case 'energy.Partial_Settings'
                cache_settings.enabled = false;
                
                [ objectiveFunction, objectiveGradient,...
                  objectiveFunction_v, objectiveGradient_v,...
                  x0_v, C_init, est_rank ] = energy.compute( Src, Tar, energy_settings );
              
                % ---
                % Fmap settings
                pair_settings = shape.PartialPairSettings(Src, Tar, singleBasisEvaluation_settings.dataset,...
                                   objectiveFunction, objectiveGradient,...
                                   objectiveFunction_v, objectiveGradient_v,...
                                   x0_v, est_rank,...
                                   'C_init', C_init);
                Fmap_settings = Fmap.Settings(energy_settings,...
                              'optimizer', 'manopt',...
                              'partial', true);
                % ---
                                               
                Fmap_info = Fmap.compute( Src, Tar, Fmap_settings, pair_settings, cache_settings );
                cache_settings.enabled = useCache;
                
                pF = Fmap_info.matches;

                % p2pMap settings
                p2pMap_settings = p2pMap.Settings(Fmap_settings,...
                                                Src.basis_settings, Tar.basis_settings,...
                                                'type' ,'KNN');
        end

        C = Fmap_info.C;
        save(sprintf('%s/%d_%s-%s_C.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'C');
        save(sprintf('%s/%d_%s-%s_pF.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'pF');
        
        if singleBasisEvaluation_settings.computeICP
            Src.BASIS_ICP = basis.compute(Src.SHAPE, basis_settings_ICP, cache_settings);
            Tar.BASIS_ICP = basis.compute(Tar.SHAPE, basis_settings_ICP, cache_settings);
            % Fmap refinement
            [C_ICP,pF_ICP] = Fmap.refine(Src, Tar, Fmap_info.C,...
                                    'type', 'ICP','parameters', refinement_parameters);
%             % ICP p2pMap computation
%             pF_ICP = p2pMap.compute(Src.BASIS, Tar.BASIS, C_ICP, p2pMap_settings );

            save(sprintf('%s/%d_%s-%s_C_ICP.mat', char(saveDir), i, shape_name_Src, shape_name_Tar),'C_ICP');
            save(sprintf('%s/%d_%s-%s_pF_ICP.mat', char(saveDir), i, shape_name_Src, shape_name_Tar),'pF_ICP');
        end;

        if singleBasisEvaluation_settings.computeZoomOut
            Src.BASIS_ZoomOut = basis.compute(Src.SHAPE, basis_settings_ZoomOut, cache_settings);
            Tar.BASIS_ZoomOut = basis.compute(Tar.SHAPE, basis_settings_ZoomOut, cache_settings);
            % Fmap refinement
            [C_ZoomOut,pF_ZoomOut] = Fmap.refine(Src, Tar, Fmap_info.C,...
                                    'type', 'ZoomOut','parameters', refinement_parameters);

            save(sprintf('%s/%d_%s-%s_C_ZoomOut.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'C_ZoomOut');
            save(sprintf('%s/%d_%s-%s_pF_ZoomOut.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'pF_ZoomOut');
        end

        if singleBasisEvaluation_settings.computeBCICP
            Src.BASIS_BCICP = basis.compute(Src.SHAPE, basis_settings_BCICP, cache_settings);
            Tar.BASIS_BCICP = basis.compute(Tar.SHAPE, basis_settings_BCICP, cache_settings);
            
            pair_settings_inverse = shape.PairSettings(Tar, Src, singleBasisEvaluation_settings.dataset, objectiveFunction, objectiveGradient);
            p2pMap_settings_inverse = p2pMap.Settings(Fmap_settings,...
                                Tar.basis_settings, Src.basis_settings,...
                                'type' ,'KNN');

            cache_settings.enabled = false;
            Fmap_info_inverse = Fmap.compute(Src, Tar, Fmap_settings, pair_settings_inverse, cache_settings);
            cache_settings.enabled = useCache;
            pF_inverse = p2pMap.compute(Tar.BASIS, Src.BASIS,...
                                Fmap_info_inverse.C, p2pMap_settings_inverse );

            [pF_BCICP, ~] = p2pMap.bcicp(Src, Tar, pF, pF_inverse, refinement_parameters.numIters_BCICP);

            C_inverse = Fmap_info_inverse.C;
            
            save(sprintf('%s/%d_%s-%s_C_inverse.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'C_inverse');
            save(sprintf('%s/%d_%s-%s_pF_BCICP.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'pF_BCICP');
            save(sprintf('%s/%d_%s-%s_pF_inverse.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'pF_inverse');
        end


        % Compute metrics
        
        % Useful data for partiality plot with the partial dataset
        switch singleBasisEvaluation_settings.dataset
            case {'Shrec16_Partial/holes', 'Shrec16_Partial/cuts'}
                result.area_Src = Src.SHAPE.sqrt_area^2;
                result.area_Tar = Tar.SHAPE.sqrt_area^2;
        end

        switch singleBasisEvaluation_settings.metricType
            case 'Euclidean'
                result.Raw = evaluation.euclidean_metrics(Src.SHAPE, pF, evalLmk_Src, evalLmk_Tar);
                if singleBasisEvaluation_settings.computeICP
                    result.ICP = evaluation.euclidean_metrics(Src.SHAPE, pF_ICP, evalLmk_Src, evalLmk_Tar);
                end
                if singleBasisEvaluation_settings.computeZoomOut
                    result.ZoomOut = evaluation.euclidean_metrics(Src.SHAPE, pF_ZoomOut, evalLmk_Src, evalLmk_Tar);
                end
                if singleBasisEvaluation_settings.computeBCICP
                    result.BCICP = evaluation.euclidean_metrics(Src.SHAPE, pF_BCICP, evalLmk_Src, evalLmk_Tar);
                end
            case 'Geodesic'
                result.Raw = evaluation.geodesic_metrics(Src.SHAPE, pF, evalLmk_Src, evalLmk_Tar);
                if singleBasisEvaluation_settings.computeICP
                    result.ICP = evaluation.geodesic_metrics(Src.SHAPE, pF_ICP, evalLmk_Src, evalLmk_Tar);
                end
                if singleBasisEvaluation_settings.computeZoomOut
                    result.ZoomOut = evaluation.geodesic_metrics(Src.SHAPE, pF_ZoomOut, evalLmk_Src, evalLmk_Tar);
                end
                if singleBasisEvaluation_settings.computeBCICP
                    result.BCICP = evaluation.geodesic_metrics(Src.SHAPE, pF_BCICP, evalLmk_Src, evalLmk_Tar);
                end
            case 'Distortion'
                result.Raw = evaluation.distortion_metrics(Src.SHAPE, Tar.SHAPE, pF);
                if singleBasisEvaluation_settings.computeICP
                    result.ICP = evaluation.distortion_metrics(Src.SHAPE, Tar.SHAPE, pF_ICP);
                end
                if singleBasisEvaluation_settings.computeZoomOut
                    result.ZoomOut = evaluation.distortion_metrics(Src.SHAPE, Tar.SHAPE, pF_ZoomOut);
                end
                if singleBasisEvaluation_settings.computeBCICP
                    result.BCICP = evaluation.distortion_metrics(Src.SHAPE, Tar.SHAPE, pF_BCICP);
                end

            case 'Map Dirichlet Energy'
                result.Raw = evaluation.map_dirichlet_energy(Src.SHAPE, Tar.SHAPE, pF);
                if singleBasisEvaluation_settings.computeICP
                    result.ICP = evaluation.map_dirichlet_energy(Src.SHAPE, Tar.SHAPE, pF_ICP);
                end
                if singleBasisEvaluation_settings.computeZoomOut
                    result.ZoomOut = evaluation.map_dirichlet_energy(Src.SHAPE, Tar.SHAPE, pF_ZoomOut);
                end
                if singleBasisEvaluation_settings.computeBCICP
                    result.BCICP = evaluation.map_dirichlet_energy(Src.SHAPE, Tar.SHAPE, pF_BCICP);
                end
        end
        
        save(sprintf('%s/%d_%s-%s_result.mat',char(saveDir), i, shape_name_Src, shape_name_Tar),'-struct','result');

    end
    
end

