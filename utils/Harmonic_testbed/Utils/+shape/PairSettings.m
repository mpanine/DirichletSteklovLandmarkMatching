classdef PairSettings < core.Settings
    %PAIRSETTINGS Settings class for a shape pair corresponding to a
    %dataset.
    
    properties
        dataset % Dataset associated with the shape pair
        name_Src % Name of the Source shape
        name_Tar % Name of the Target shape
        basis_settings_Src % Basis settings used to compute the basis on the Source
        basis_settings_Tar % Basis settings used to compute the basis on the Target
        name % Name of the shape pair
        
        constantFunction % Value of the first column of a Fmap computed with
        % MHB for Source and Target shapes (used for the default projection function)
        C_init % Initial value for which to begin the Fmap estimation
        objectiveFunction % Objective function used in the Fmap estimation
        objectiveGradient % Gradient of the objective function used in the Fmap estimation
        projectionFunction % Projection function used in the Fmap estimation
    end
    
    methods
        function settings = PairSettings(Src, Tar, dataset, objectiveFunction, objectiveGradient, varargin)
            assert(isa(objectiveFunction, 'function_handle'),...
                      sprintf('objectiveFunction has to be a function_handle object got %s object instead.',class(objectiveFunction)));
            assert(isa(objectiveGradient, 'function_handle'),...
                      sprintf('objectiveGradient has to be a function_handle object got %s object instead.',class(objectiveGradient)));
            settings.dataset = dataset;
            settings.name_Src = Src.SHAPE.name;
            settings.name_Tar = Tar.SHAPE.name;
            settings.basis_settings_Src = Src.basis_settings;
            settings.basis_settings_Tar = Tar.basis_settings;
            
            constantFunction_default = zeros(Tar.basis_settings.numBasisFun,1);
            constantFunction_default(1) = sign(Src.BASIS.basis(1,1)*Tar.BASIS.basis(1,1))*sqrt(sum(Tar.SHAPE.area)/sum(Src.SHAPE.area));
            
            settings.constantFunction = constantFunction_default;
        
            settings.objectiveFunction = objectiveFunction;
            settings.objectiveGradient = objectiveGradient;
            
            settings.name = sprintf('%s %s',...
                                settings.name_Src, settings.name_Tar);
                            
            settings.objectiveFunction = objectiveFunction;
            
            
            C_init_default = zeros(settings.basis_settings_Tar.numBasisFun*...
                           settings.basis_settings_Src.numBasisFun, 1);
            C_init_default(1) = settings.constantFunction(1);
            
            funProj_default = @(F) [settings.constantFunction; ...
                F(settings.basis_settings_Tar.numBasisFun+1:end)];
            
            p = inputParser;
            addParameter(p,'C_init', C_init_default, @isnumeric);
            addParameter(p,'projectionFunction', funProj_default, @isnumeric);
            parse(p,varargin{:});
            inputs = p.Results;
            
            
            settings.projectionFunction = inputs.projectionFunction;
            settings.C_init = inputs.C_init;
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
        
        function self = set.dataset(self, value)
            assert(ischar(value), 'dataset has to be a string.');
            self.dataset = value;
            self = self.update();
        end
        
        function self = set.name(self, value)
            assert(ischar(value), 'name has to be a string.');
            self.name = value;
            self = self.update();
        end
        
        function self = set.name_Src(self, value)
            assert(ischar(value), 'name_Src has to be a string.');
            self.name_Src = value;
            self = self.update();
        end
        
        function self = set.name_Tar(self, value)
            assert(ischar(value), 'nameTar has to be a string.');
            self.name_Tar = value;
            self = self.update();
        end
        
        function self = set.constantFunction(self, value)
            assert(isnumeric(value), 'constantFunction has to be a numerical array.');
            self.constantFunction = value;
            self = self.update();
        end
        
        function self = set.objectiveFunction(self, value)
            assert(isa(value, 'function_handle'), 'objectiveFunction has to be a function handle.');
            self.objectiveFunction = value;
            self = self.update();
        end
        
        function self = set.objectiveGradient(self, value)
            assert(isa(value, 'function_handle'), 'objectiveGradient has to be a function handle.');
            self.objectiveGradient = value;
            self = self.update();
        end
        
        function self = set.projectionFunction(self, value)
            assert(isa(value, 'function_handle'), 'projectionFunction has to be a function handle.');
            self.projectionFunction = value;
            self = self.update();
        end
        
        function self = set.C_init(self, value)
            assert(isnumeric(value), 'C_init has to be a numerical array.');
            self.C_init = value;
            self = self.update();
        end
    end
    
end

