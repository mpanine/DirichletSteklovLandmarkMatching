classdef LMH_Settings < basis.Settings
    %LMH_SETTINGS Setting for the LMH basis.
    
    properties
    end
    
    methods
        function settings = LMH_Settings(numBasisFun, k_local, mu, eta)
            if nargin == 0
                numBasisFun = 50;
                k_local = 25;
                mu = 1e-2;
                eta = 1e5;
                warning('Using default values for basis.LMH_Settings initialization.\nDefault values:\n\t numBasisFun=%d\n\t k_local=%d\n\t mu=%d\n\t eta=%d',numBasisFun,k_local,mu,eta);
            end
            assert(isfinite(k_local) & k_local==floor(k_local),...
                'k_local has to be an integer.');
            assert(isnumeric(mu), 'mu has to be a numerical value.');
            assert(isnumeric(eta), 'eta has to be a numerical value.');
            
            parameters.k_local = k_local; % The number of localized manifold harmonic functions to compute.
            % k_local <= numBasisFun
            % if k_local < numBasisFun, the first (numBasisFun - k_local) MHB functions will complete the basis.
            parameters.mu = mu; % Controles the locality of the basis.
            parameters.eta = eta; % Controles the orthogonality of the basis.
            
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'LMH';
            
            settings = settings.update();
            
        end
    end
    
end

