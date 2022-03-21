classdef CMM_Settings < basis.Settings
    %CMM_SETTINGS Settings for the CMM basis
    
    properties
    end
    
    methods
        function settings = CMM_Settings(numBasisFun, mu, accelerated)
            if nargin == 0
                numBasisFun = 50;
                mu = 1e-2;
                accelerated = true;
                warning('Using default values for basis.CMM_Settings initialization.\nDefault values:\n\t numBasisFun=%d\n\t mu=%d\n\t accelerated=%d',numBasisFun,mu,accelerated);
            end
            assert(isnumeric(mu), 'mu has to be a numerical value.');
            
            if nargin < 3
                accelerated = true;
            end;
            
            parameters = struct;
            parameters.mu = mu; % Controles the locality of the basis.
            parameters.accelerated = accelerated;
            
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'CMM';
            
            settings = settings.update();
        end
    end
    
end

