classdef Kernel_Settings < descriptor.Settings
    %KERNEL_SETTINGS Settings class for Wave/Heat kernel descriptors
    
    properties
%         Fields of Kernel_Settings.parameters:
%         numSkip % One out of numSkip descriptors are taken to reduce the number of descriptors.
%         numTimes % The number of times used to compute the descriptors. Is equal to the number of descriptors computed before downsampling.
%         numEigs % The number of eigenpairs used to compute the descriptors.
%         regions % Not used for now.
%         regionsWeights % Not used for now.
%         functionType % The type of descriptor function to compute:
%          - 'wave': wave kernel signature descriptors.
%          - 'heat': heat kernel signature descriptors.
%         computeSignatureFun % If true, adds the signature function associated to the functionType at the begining of the descriptors array (before downsampling).
        
        haveRegions % Is true when regions ~= []
        haveWeights % Is true when regionWeights ~= []
    end
    
    methods
        function settings = Kernel_Settings(shape_settings, varargin)
            
            p = inputParser;
            addParameter(p,'landmarks', [], @isnumeric);
            addParameter(p,'numSkip', 40, @isnumeric);
            addParameter(p,'numTimes', 200, @isnumeric);
            addParameter(p,'numEigs',...
                min(shape_settings.MHB_settings.numBasisFun, 100), @isnumeric);
            addParameter(p,'regions',{},@iscell);
            addParameter(p,'regionsWeights', [], @isnumeric);
            addParameter(p,'functionType', 'heat', @ischar);
            addParameter(p,'computeSignatureFun', true, @islogical);
            parse(p,varargin{:});
            inputs = p.Results;
            
            parameters.numSkip = inputs.numSkip;
            parameters.numTimes = inputs.numTimes;
            parameters.numEigs = inputs.numEigs;
            parameters.regions = inputs.regions;
            parameters.regionsWeights = inputs.regionsWeights;
            parameters.functionType = inputs.functionType;
            parameters.computeSignatureFun = inputs.computeSignatureFun;
            
            settings@descriptor.Settings(inputs.landmarks, parameters);
            
            settings.name = 'Kernel';

            if isempty(settings.parameters.regions)
                settings.haveRegions = false;
            else
                settings.haveRegions = true;
            end

            if isempty(settings.parameters.regionsWeights)
                settings.haveWeights = false;
            else
                settings.haveWeights = true;
            end
            
            settings = settings.update();
        end
        
        function self = set.haveRegions(self, value)
            assert(islogical(value), 'haveRegions has to be a logical.');
            self.haveRegions = value;
            self = self.update();
        end
        
        function self = set.haveWeights(self, value)
            assert(islogical(value), 'haveWeights has to be a logical.');
            self.haveWeights = value;
            self = self.update();
        end
        
%         function self = set.numSkip(self, value)
%             assert(isfinite(value) & value==floor(value), 'numSkip has to be an integer.');
%             self.numSkip = value;
%             self = self.update();
%         end
%         
%         function self = set.numTimes(self, value)
%             assert(isfinite(value) & value==floor(value), 'numTimes has to be an integer.');
%             self.numTimes = value;
%             self = self.update();
%         end
%         
%         function self = set.numEigs(self, value)
%             assert(isfinite(value) & value==floor(value), 'numEigs has to be an integer.');
%             self.numEigs = value;
%             self = self.update();
%         end
%         
%         function self = set.regions(self, value)
%             assert(iscell(value), 'regions has to be a cell array.');
%             self.regions = value;
%             self = self.update();
%         end
%         
%         function self = set.regionsWeights(self, value)
%             assert(isnumeric(value), 'regionsWeights has to be a numerical array.');
%             self.regionsWeights = value;
%             self = self.update();
%         end
%         
%         function self = set.functionType(self, value)
%             assert(ischar(value), 'functionType has to be a string.');
%             self.functionType = value;
%             self = self.update();
%         end
%         
%         function self = set.computeSignatureFun(self, value)
%             assert(islogical(value), 'computeSignatureFun has to be a logical.');
%             self.computeSignatureFun = value;
%             self = self.update();
%         end
    end
    
end
