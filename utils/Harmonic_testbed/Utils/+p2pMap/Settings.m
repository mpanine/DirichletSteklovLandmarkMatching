classdef Settings < core.Settings
    %SETTINGS Settings for the point to point map computation
    
    properties
        type % Type of p2p map computation
        useBCICP % whether to use or not BCICP
        Fmap_settings % Fmap.Settings used to compute the Fmap 
                      % to be converted to a p2p map
        basis_settings_Src % basis.Settings object for the Source shape 
                           % to be used to compute the p2p map
        basis_settings_Tar % basis.Settings object for the Target shape 
                           % to be used to compute the p2p map
    end
    
    methods
        function settings = Settings(Fmap_settings, basis_settings_Src, basis_settings_Tar, varargin)
            p = inputParser;
            addParameter(p,'type','KNN',@ischar);
            addParameter(p,'useBCICP',false,@islogical);
            parse(p,varargin{:});
            inputs = p.Results;
            
            settings.type = inputs.type;
            settings.useBCICP = inputs.useBCICP;
            settings.Fmap_settings = Fmap_settings;
            settings.basis_settings_Src = basis_settings_Src;
            settings.basis_settings_Tar = basis_settings_Tar;
        end
        
        function self = set.type(self, value)
            assert(ischar(value), 'type has to be a string.');
            self.type = value;
            self = self.update();
        end
        
        function self = set.Fmap_settings(self, value)
            assert(isa(value, 'Fmap.Settings'),...
                      sprintf('Fmap_settings has to be a Fmap.Settings object got %s object instead.',class(value)));
            self.Fmap_settings = value;
            self = self.update();
        end
        
        function self = set.basis_settings_Src(self, value)
            assert(isa(value, 'basis.Settings'),...
                      sprintf('basis_settings_Src has to be a basis.Settings object got %s object instead.',class(value)));
            self.basis_settings_Src = value;
            self = self.update();
        end
        
        function self = set.basis_settings_Tar(self, value)
            assert(isa(value, 'basis.Settings'),...
                      sprintf('basis_settings_Tar has to be a basis.Settings object got %s object instead.',class(value)));
            self.basis_settings_Tar = value;
            self = self.update();
        end
    end
    
end

