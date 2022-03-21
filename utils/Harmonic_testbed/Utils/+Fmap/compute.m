function [ Fmap_info ] = compute( Src, Tar, Fmap_settings, pair_settings, cache_settings, options )
%COMPUTE Computes the Fmap associated to the Source and Target given the
%Fmap_settings.
%   Inputs:
%       - pair_settings: shape.PairSettings object
%       - Fmap_settings: Fmap.Settings object
%       - cache_settings: cache.Settings object
%       - options: struct of optimizer options

    if nargin < 6
        options.verbose = 1;
        options.maxIter = 1e6;
    end
    
    partial_str = '';
    if Fmap_settings.partial
        partial_str = 'partial ';
    end
    fprintf('Computing %sFmap...',partial_str); tic;

    assert(isa(Fmap_settings, 'Fmap.Settings'),...
        sprintf('Fmap_settings has to be a Fmap.Settings object got %s object instead.',...
                                                    class(Fmap_settings)));
    assert(isa(pair_settings, 'shape.PairSettings'),...
        sprintf('pair_settings has to be a shape.PairSettings object got %s object instead.',...
                                                    class(pair_settings)));
    assert(isa(cache_settings, 'cache.Settings'),...
        sprintf('cache_settings has to be a cache.Settings object got %s object instead.',...
                                                    class(cache_settings)));

    cached_Fmap = cache.load(cache_settings, pair_settings.name, Fmap_settings);
    
    if length(fieldnames(cached_Fmap)) == 0
        numBasisFun_Src = pair_settings.basis_settings_Src.numBasisFun;
        numBasisFun_Tar = pair_settings.basis_settings_Tar.numBasisFun;
        
        if Fmap_settings.partial
            Fmap_info = Fmap.compute_partial(Src, Tar, Fmap_settings, pair_settings, cache_settings, options);
        else
            switch Fmap_settings.optimizer
                case 'manopt'
                    manifold = euclideanfactory(numBasisFun_Tar,...
                                                numBasisFun_Src);
                    problem = {};

                    problem.M = manifold;

                    problem.cost = pair_settings.objectiveFunction;
                    problem.egrad = pair_settings.objectiveGradient;

                    options.verbosity = 0;
                    [C_vec, final_cost, info, ~] = conjugategradient(problem,...
                                                pair_settings.C_init, options);
                    C = reshape(C_vec,[numBasisFun_Tar,...
                                       numBasisFun_Src]);
                    fprintf('Manopt optimization, iterations: %d, cost: %f\n',...
                            info(end).iter, final_cost);
                otherwise
                    funObj = @(C) deal(pair_settings.objectiveFunction(C),...
                                       pair_settings.objectiveGradient(C));
                    C = reshape(minConf_PQN(funObj,...
                                            pair_settings.C_init,...
                                            pair_settings.projectionFunction,...
                                            options),...
                                            [numBasisFun_Tar,...
                                             numBasisFun_Src]);
            end
            Fmap_info.settings = Fmap_settings;
            Fmap_info.pair_settings = pair_settings;
            Fmap_info.C = C;
        end
                             
        cache.save(Fmap_info, pair_settings.name, cache_settings);
    else
        Fmap_info = cached_Fmap.object;
    end
    
    t =toc; fprintf('done:%.4fs\n',t);
end

