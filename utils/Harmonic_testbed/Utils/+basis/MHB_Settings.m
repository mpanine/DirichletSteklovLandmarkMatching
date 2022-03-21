classdef MHB_Settings < basis.Settings
    %MHB_SETTINGS Settings for the MHB basis
    
    properties
    end
    
    methods
        function settings = MHB_Settings(numBasisFun)
            if nargin == 0
                numBasisFun = 50;
                warning('Using default values for basis.MHB_Settings initialization.\nDefault values:\n\t numBasisFun=%d',numBasisFun);
            end
            parameters = struct;
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'MHB';
            
            settings = settings.update();
        end
    end
    
end

