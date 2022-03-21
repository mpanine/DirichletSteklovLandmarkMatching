classdef LDB_Settings < basis.Settings
    %LDB_SETTINGS Settings for the LDB basis
    
    properties
    end
    
    methods
        function settings = LDB_Settings(mu, numFunPerMu)
            if nargin == 0
                mu = [1e-3,1e-2,1e-1,1];
                numFunPerMu = 12;
                warning('Using default values for basis.LDB_Settings initialization.\nDefault values:\n\t mu=[%d,%d,%d,%d]\n\t numFunPerMu=%d',mu,numFunPerMu);
            end
            
            assert(isfinite(numFunPerMu) & numFunPerMu==floor(numFunPerMu),...
                'numFunPerMu has to be an integer.');
            assert(isnumeric(mu), 'mu has to be a numerical list.');
            
            parameters = struct;
            parameters.numFunPerMu = numFunPerMu; % Number of basis functions computed for each value of mu.
            parameters.mu = mu; % List of diffusion times used to compute the basis.
            
            numBasisFun = numFunPerMu*length(mu);
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'LDB';
            
            settings = settings.update();
        end
    end
    
end

