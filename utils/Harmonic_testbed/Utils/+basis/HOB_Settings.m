classdef HOB_Settings < basis.Settings
    %HOB_SETTINGS Settings for the MHB basis
    
    properties
    end
    
    methods
        function settings = HOB_Settings(numBasisFun)
            if nargin == 0
                numBasisFun = 50;
                warning('Using default values for basis.HOB_Settings initialization.\nDefault values:\n\t numBasisFun=%d',numBasisFun);
            end
            parameters = struct;
%             parameters.depth = depth; % The number of neighboring vertices to use when computing.
            
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'HOB';
            
            settings = settings.update();
        end
    end
    
end

