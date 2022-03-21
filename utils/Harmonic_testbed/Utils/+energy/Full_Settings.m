classdef Full_Settings < energy.Settings
    %FULL_SETTINGS Settings for the energy computation in the case of full Fmaps.
    
    properties
        computeLBCommutativity
        computeDescriptorPreservation
        computeDescriptorCommutativity
        computeOrientationTerm
        
        useResolvent
        useOperatorNormalization;
        useLegacyNormalization;
        
        weight_LBCommutativity
        weight_descriptorPreservation
        weight_descriptorCommutativity
        weight_orientationTerm
    end
    
    methods
        function settings = Full_Settings(basis_settings_Src, ...
                                     basis_settings_Tar, ...
                                     varargin)
                                 
             settings@energy.Settings(basis_settings_Src,...
                                      basis_settings_Tar);
            
            p = inputParser;
            addParameter(p,'weight_LBCommutativity', 1e-3, @isnumeric);
            addParameter(p,'weight_descriptorCommutativity', 1, @isnumeric);
            addParameter(p,'weight_descriptorPreservation', 1e-1, @isnumeric);
            addParameter(p,'weight_orientationTerm', 0, @isnumeric);
            addParameter(p, 'useResolvent', false, @islogical);
            addParameter(p, 'useOperatorNormalization', false, @islogical);
            addParameter(p, 'useLegacyNormalization', false, @islogical);
            parse(p,varargin{:});
            inputs = p.Results;
            
            settings.weight_LBCommutativity = inputs.weight_LBCommutativity;
            settings.weight_descriptorPreservation = inputs.weight_descriptorPreservation;
            settings.weight_descriptorCommutativity = inputs.weight_descriptorCommutativity;
            settings.weight_orientationTerm = inputs.weight_orientationTerm;
            
            settings.useResolvent = inputs.useResolvent;
            settings.useOperatorNormalization = inputs.useOperatorNormalization;
            settings.useLegacyNormalization = inputs.useLegacyNormalization;

            if settings.weight_LBCommutativity == 0
                settings.computeLBCommutativity = false;
            else
                settings.computeLBCommutativity = true;
            end

            if settings.weight_descriptorPreservation == 0
                settings.computeDescriptorPreservation = false;
            else
                settings.computeDescriptorPreservation = true;
            end

            if settings.weight_descriptorCommutativity == 0
                settings.computeDescriptorCommutativity = false;
            else
                settings.computeDescriptorCommutativity = true;
            end

            if settings.weight_orientationTerm == 0
                settings.computeOrientationTerm = false;
            else
                settings.computeOrientationTerm = true;
            end
            
        end
        
        function self = set.weight_LBCommutativity(self, value)
            assert(isnumeric(value), 'weight_LBCommutativity has to be a numerical value.');
            self.weight_LBCommutativity = value;
            self = self.update();
        end
        
        function self = set.weight_descriptorPreservation(self, value)
            assert(isnumeric(value), 'weight_descriptorPreservation has to be a numerical value.');
            self.weight_descriptorPreservation = value;
            self = self.update();
        end
        
        function self = set.weight_descriptorCommutativity(self, value)
            assert(isnumeric(value), 'weight_descriptorCommutativity has to be a numerical value.');
            self.weight_descriptorCommutativity = value;
            self = self.update();
        end

        function self = set.weight_orientationTerm(self, value)
            assert(isnumeric(value), 'weight_orientationTerm has to be a numerical value.');
            self.weight_orientationTerm = value;
            self = self.update();
        end
        
        function self = set.useResolvent(self, value)
            assert(islogical(value), 'useResolvent has to be a logical value.');
            self.useResolvent = value;
            self = self.update();
        end
        
        function self = set.useOperatorNormalization(self, value)
            assert(islogical(value), 'useOperatorNormalization has to be a logical value.');
            self.useOperatorNormalization = value;
            self = self.update();
        end
        
        function self = set.useLegacyNormalization(self, value)
            assert(islogical(value), 'useLegacyNormalization has to be a logical value.');
            self.useLegacyNormalization = value;
            self = self.update();
        end
    end
    
end

