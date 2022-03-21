classdef MWB_Settings < basis.Settings
    %MWB_SETTINGS Settings for the MHB basis
    
    properties
    end
    
    methods
        function settings = MWB_Settings(numWaveletsPerScale,timeScales,type,orthonormalization)
            if nargin == 0
                numLevels = 10;
                log_ts = linspace(log(1e-3), log(1e0), numLevels);
                timeScales = exp(log_ts);
                numWaveletsPerScale = 1;
                orthonormalization = false;
                type = 'mexican hat'
                warning('Using default values for basis.MWB_Settings initialization.\nDefault values:\n\t numWaveletsPerScale=%d\n\t timeScales=%2f\n\t type=%s\n\t orthonormalization=%d',numWaveletsPerScale,timeScales,type,orthonormalization);
            end
            
            parameters = struct;
            parameters.timeScales = timeScales;
            parameters.numWaveletsPerScale = numWaveletsPerScale;
            parameters.type = type;
            parameters.orthonormalization = orthonormalization;
            numBasisFun = length(parameters.timeScales)*parameters.numWaveletsPerScale;
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'MWB';
            
            settings = settings.update();
        end
    end
    
end

