classdef Partial_Settings < energy.Settings
    %PARTIAL_SETTINGS Settings contains the additional settings used for 
    %partial Fmaps.
    
    properties
        weight_slantedDiagonalMask % slanted-diagonal mask (mu1 in original code)
        weight_subOrthogonality % sub-orthogonality (mu2 in original code)
        weight_area % area (mu3 in original code)
        weight_regularity % regularity (mu4 in original code)
        
        tv_sigma
        tv_mean
    end
    
    methods
        function settings = Partial_Settings( basis_settings_Src, basis_settings_Tar, varargin )
            
            settings@energy.Settings(basis_settings_Src,...
                                      basis_settings_Tar);
                                  
            p = inputParser;
            addParameter(p,'weight_slantedDiagonalMask', 1, @isnumeric);
            addParameter(p,'weight_subOrthogonality', 1e3, @isnumeric);
            addParameter(p,'weight_area', 1, @isnumeric);
            addParameter(p,'weight_regularity', 1e2, @isnumeric);
            
            % TOSCA parameters by default
            addParameter(p,'tv_sigma', 0.2, @isnumeric);
            addParameter(p,'tv_mean', 0.5, @isnumeric);
            parse(p,varargin{:});
            inputs = p.Results;
            
            settings.weight_slantedDiagonalMask = inputs.weight_slantedDiagonalMask;
            settings.weight_subOrthogonality = inputs.weight_subOrthogonality;
            settings.weight_area = inputs.weight_area;
            settings.weight_regularity = inputs.weight_regularity;
            settings.tv_sigma = inputs.tv_sigma;
            settings.tv_mean = inputs.tv_mean;
        end
        
        function self = set.weight_slantedDiagonalMask(self, value)
            assert(isnumeric(value), 'weight_slantedDiagonalMask has to be numerical.');
            self.weight_slantedDiagonalMask = value;
            self = self.update();
        end
        
        function self = set.weight_subOrthogonality(self, value)
            assert(isnumeric(value), 'weight_subOrthogonality has to be numerical.');
            self.weight_subOrthogonality = value;
            self = self.update();
        end
        
        function self = set.weight_area(self, value)
            assert(isnumeric(value), 'weight_area has to be numerical.');
            self.weight_area = value;
            self = self.update();
        end
        
        function self = set.weight_regularity(self, value)
            assert(isnumeric(value), 'weight_regularity has to be numerical.');
            self.weight_regularity = value;
            self = self.update();
        end
        
        function self = set.tv_sigma(self, value)
            assert(isnumeric(value), 'tv_sigma has to be numerical.');
            self.tv_sigma = value;
            self = self.update();
        end
        
        function self = set.tv_mean(self, value)
            assert(isnumeric(value), 'tv_mean has to be numerical.');
            self.tv_mean = value;
            self = self.update();
        end
        
    end
    
end

