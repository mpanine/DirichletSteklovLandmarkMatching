classdef Settings < core.Settings
    %SETTINGS Settings for the energy computation base class.
    
    properties
        
        basis_settings_Src
        basis_settings_Tar
    end
    
    methods
        function settings = Settings(basis_settings_Src, ...
                                     basis_settings_Tar)
                                 
            assert(isa(basis_settings_Src, 'basis.Settings'),...
                   sprintf('basis_settings_Src has to be a basis.Settings object got %s object instead.',class(basis_settings_Src)));
            assert(isa(basis_settings_Tar, 'basis.Settings'),...
                   sprintf('basis_settings_Tar has to be a basis.Settings object got %s object instead.',class(basis_settings_Tar)));
            
            settings.basis_settings_Src = basis_settings_Src;
            settings.basis_settings_Tar = basis_settings_Tar;
            
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

