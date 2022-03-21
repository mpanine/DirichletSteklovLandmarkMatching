classdef Settings < core.Settings
    %SETTINGS Superclass for all basis settings classes
    
    properties
        numBasisFun % The number of basis functions in the basis.
        parameters % The parameters that are required to compute the basis.
        name % The name of the basis. It serves as an identifier.
    end
    
    methods
        function settings = Settings(numBasisFun, parameters)
            if nargin < 2
                parameters = struct;
            end
            settings.numBasisFun = numBasisFun;
            settings.parameters = parameters; % Struct containing all additional
            % computing parameters of the basis. Empty if no parameters.
            
            settings.name = 'None';
        end
        
        function self = set.numBasisFun(self, value)
            assert(isfinite(value) & value==floor(value), 'numBasisFun has to be an integer.');
            self.numBasisFun = value;
            self = self.update();
        end
        
        function self = set.name(self, value)
            assert(ischar(value), 'name has to be a string.');
            self.name = value;
            self = self.update();
        end
        
        function self = set.parameters(self, value)
            assert(isstruct(value), 'parameters has to be a struct.');
            self.parameters = value;
            self = self.update();
        end
    end
    
end

