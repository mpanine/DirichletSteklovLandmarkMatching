classdef Shot_Settings < descriptor.Settings
    %SHOT_SETTINGS Settings for Shot descriptors.
    
    properties
%         Fields of Kernel_Settings.parameters:
%             - numBins % The number of bins to compute
%             - radius % The radius of the Shot descriptors
%             - numNeighbors % The number of neighbors used in the
%             computation
%             - numSkip % The number of times used to compute the descriptors. Is equal to the number of descriptors computed before downsampling.
    end
    
    methods
        function settings = Shot_Settings(varargin)
            
            p = inputParser;
            addParameter(p,'landmarks', [], @isnumeric);
            addParameter(p,'numBins', 9, @isnumeric);
            addParameter(p,'radius', 0.01, @isnumeric);
            addParameter(p,'numNeighbors', 3, @isnumeric);
            addParameter(p,'numSkip', 1, @isnumeric);
            parse(p,varargin{:});
            inputs = p.Results;
            
            parameters.numSkip = inputs.numSkip;
            parameters.numBins = inputs.numBins;
            parameters.numNeighbors = inputs.numNeighbors;
            parameters.radius = inputs.radius;
            
            settings@descriptor.Settings(inputs.landmarks, parameters);
            
            settings.name = 'Shot';
            
            settings = settings.update();
        end
    end
    
end

