% original file: test_commute_faust.m
% paper: Informative Descriptor Preservation via Commutativity for Shape Matching
function [ DESCRIPTORS ] = compute(S, descriptor_settings, cache_settings)
    A = S.A;
    
    fct = [];

    DESCRIPTORS.settings = descriptor_settings;

    cached_descriptors = cache.load(cache_settings, S.name, descriptor_settings);

    if length(fieldnames(cached_descriptors)) == 0
        switch class(descriptor_settings)
            case 'descriptor.Kernel_Settings'
                numEigs = descriptor_settings.parameters.numEigs;
                numTimes = descriptor_settings.parameters.numTimes;
                
                numMapFun = length(descriptor_settings.landmarks);
                fprintf('Computing %d %s kernel descriptor functions...',...
                    numTimes/descriptor_settings.parameters.numSkip,...
                    descriptor_settings.parameters.functionType); tic;

                Basis = S.laplacianBasis.basis(:,1:numEigs);
                Ev = S.laplacianBasis.eigenvalues(1:numEigs);

                % only skip the WKS descriptors but keep all the lmk/region descriptors
                if descriptor_settings.parameters.computeSignatureFun
                    switch descriptor_settings.parameters.functionType
                        case 'wave'
                            fct = [fct, waveKernelSignature(Basis, Ev, A, numTimes)];
                        case 'heat'
                            fct = [fct, heatKernelSignature(Basis, Ev, A, numTimes)];
                    end;
                end;

                % descriptors based on landmarks;
                if descriptor_settings.haveLmks
                    switch descriptor_settings.parameters.functionType
                        case 'wave'
                            fct = [fct, waveKernelMap(Basis, Ev, A, numTimes,...
                                                      descriptor_settings.landmarks)];
                        case 'heat'
                            fct = [fct, heatKernelMap(Basis, Ev, A, numTimes,...
                                                      descriptor_settings.landmarks)];
                    end;
                end
                fct = fct(:,1:descriptor_settings.parameters.numSkip:end);

            case 'descriptor.Shot_Settings'
                fprintf('Computing %d Shot descriptor functions...',...
                    descriptor_settings.parameters.numBins/descriptor_settings.parameters.numSkip); tic;
                
                if S.nv ~= length(descriptor_settings.landmarks)
                    throw(MException('DESCRIPTORCOMPUTE:notEnoughLandmarks',...
                                 'To compute shot descriptors, the descriptors_settings.landmarks must be equal to 1:shape.nv.'));
                end
                fct = [fct, calc_shot(S.surface.VERT', S.surface.TRIV',...
                                      descriptor_settings.landmarks,...
                                      descriptor_settings.parameters.numBins,...
                                      descriptor_settings.parameters.radius,...
                                      descriptor_settings.parameters.numNeighbors)'];
                fct = fct(:,1:descriptor_settings.parameters.numSkip:end);
            otherwise
                throw(MException('DESCRIPTORCOMPUTE:settingsUnknown',...
                                 sprintf(...
                                 'The descriptor.Settings %s is unknown.',...
                                 class(descriptor_settings))));
        end

        DESCRIPTORS.fct = fct;

        cache.save(DESCRIPTORS, S.name, cache_settings);
    else
        DESCRIPTORS = cached_descriptors.object;
    end
    t =toc; fprintf('done:%.4fs\n',t);
end
