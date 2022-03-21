classdef Settings < core.Settings
    %SETTINGS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name % The name of the descriptor
        parameters % The descriptor-specific parameters
        landmarks % The list of landmark indices where to compute the descriptors.
        haveLmks % Is true when landmarks ~= []
    end
    
    methods
        function settings = Settings(landmarks, parameters)
               
            if nargin < 2
                parameters = struct;
            end
            settings.landmarks = landmarks;
            settings.parameters = parameters; % Struct containing all additional
            % computing parameters of the descriptors. Empty if no parameters.
            
            settings.name = 'None';
            
            settings.haveLmks = false;
            if ~isempty(settings.landmarks)
                settings.haveLmks = true;
            end
            
            settings = settings.update();
        end
        
        function self = set.landmarks(self, value)
            assert(isnumeric(value), 'landmarks has to be a list of landmark indices.');
            self.landmarks = value;
            
            self = self.update_haveLmks();
            self = self.update();
        end
        
        function self = set.haveLmks(self, value)
            assert(islogical(value), 'haveLmks has to be logical.');
            self.haveLmks = value;
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
        
        function self = update_haveLmks(self)
            self.haveLmks = false;
            if ~isempty(self.landmarks)
                self.haveLmks = true;
            end
        end
    end
    
end

