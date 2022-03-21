function [ BASIS ] = compute( S, basis_settings, cache_settings, varargin )
%COMPUTE Computes the basis based on its name.
%   Input:
%   - S: the shape on which to compute the basis functions
%   - basis_settings: the basis.Settings object
%   - cache_settings: the cache.Settings object

    p = inputParser;
    addOptional(p,'harmonic_target', S);
    parse(p,varargin{:});
    harmonic_target = p.Results.harmonic_target;

    assert(isa(basis_settings, 'basis.Settings'), ...
    sprintf('basis_settings has to be of type basis.Settings. Got %s instead.', ...
    class(basis_settings)));
    assert(isa(cache_settings, 'cache.Settings'), ...
    sprintf('cache_settings has to be of type cache.Settings. Got %s instead.', ...
    class(cache_settings)));

    fprintf('Computing %d %s functions...',basis_settings.numBasisFun, basis_settings.name); tic;

    cache_name = S.name;
    if isa(basis_settings,'basis.HMB_Settings')
        cache_name = [S.name '-' harmonic_target.name];
    end
    cached_basis = cache.load(cache_settings, cache_name, basis_settings);

    if length(fieldnames(cached_basis)) == 0
        
        compute_basis = str2func(sprintf('basis.compute_%s',basis_settings.name));
        if isa(basis_settings,'basis.HMB_Settings')
            BASIS = compute_basis(S, harmonic_target, basis_settings);
        else
            BASIS = compute_basis(S, basis_settings);
        end
        BASIS.settings = basis_settings;
        
        cache.save(BASIS, cache_name, cache_settings);
    else
        BASIS = cached_basis.object;
    end
    t =toc; fprintf('done:%.4fs\n',t);
end
