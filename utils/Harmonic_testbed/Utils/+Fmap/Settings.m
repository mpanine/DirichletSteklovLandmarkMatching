classdef Settings < core.Settings
    %SETTINGS Settings for Fmap computation
    
    properties
        energy_settings % energy.Settings for the computation
        partial % To enable the partial functional map computation
        optimizer % The optimizer used when evaluating the Fmap
    end
    
    methods
        function settings = Settings(energy_settings,...
                                     varargin)
            assert(isa(energy_settings, 'energy.Settings')||isa(energy_settings, 'energy.PartialSettings'),...
                      sprintf('energy_settings has to be an energy.Settings or energy.PartialSettings object got %s object instead.',class(energy_settings)));
            p = inputParser;
            
            addParameter(p,'partial', false, @islogical);
            addParameter(p,'optimizer', 'standard', @ischar);
            parse(p,varargin{:});
            inputs = p.Results;
            
            settings.energy_settings = energy_settings;
            settings.partial = inputs.partial;
            settings.optimizer = inputs.optimizer;
        end
        
        
        function self = set.energy_settings(self, value)
            assert(isa(value, 'energy.Settings')||isa(value, 'energy.PartialSettings'),...
                      sprintf('energy_settings has to be an energy.Settings or energy.PartialSettings object got %s object instead.',class(value)));
            self.energy_settings = value;
            self = self.update();
        end
        
        function self = set.partial(self, value)
            assert(islogical(value), 'partial has to be a logical.');
            self.partial = value;
            self = self.update();
        end
        
        function self = set.optimizer(self, value)
            assert(ischar(value), 'optimizer has to be a char.');
            self.optimizer = value;
            self = self.update();
        end
    end
    
end

