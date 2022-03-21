classdef PartialPairSettings < shape.PairSettings
    %PARTIALPAIRSETTINGS Settings class for a shape pair corresponding to a
    %dataset in the partial setup.
    
    properties
        ref_BASIS_Src % basis used to refine the Fmap (ICP) on the Source shape
        ref_BASIS_Tar % basis used to refine the Fmap (ICP) on the Target shape
        objectiveFunction_v % Objective function used in the Fmap estimation for v
        objectiveGradient_v % Gradient of the objective function used in the Fmap estimation for v
        x0_v % x0 value for v optimization
        est_rank % estimated rank of the Fmap.
    
        max_iter
        icp_max_iter
        refine_fps
        fps_variance
    end
    
    methods
        function settings = PartialPairSettings(Src, Tar, dataset, objectiveFunction, objectiveGradient, objectiveFunction_v, objectiveGradient_v, x0_v, est_rank, varargin)
            settings = settings@shape.PairSettings(Src, Tar, dataset, objectiveFunction, objectiveGradient, varargin{:});
            
            assert(isa(objectiveFunction_v, 'function_handle'),...
                      sprintf('objectiveFunction_v has to be a function_handle object got %s object instead.',class(objectiveFunction_v)));
            assert(isa(objectiveGradient_v, 'function_handle'),...
                      sprintf('objectiveGradient_v has to be a function_handle object got %s object instead.',class(objectiveGradient_v)));
            assert(isa(x0_v, 'function_handle'),...
                      sprintf('x0_v has to be a function_handle object got %s object instead.',class(x0_v)));
                  
            settings.ref_BASIS_Src = Src.BASIS;
            settings.ref_BASIS_Tar = Tar.BASIS;
            
            settings.objectiveFunction_v = objectiveFunction_v;
            settings.objectiveGradient_v = objectiveGradient_v;
            settings.x0_v = x0_v;
            settings.est_rank = est_rank;
            
            p = inputParser;
            addParameter(p,'max_iter', 6, @isnumeric);
            addParameter(p,'icp_max_iter', 30, @isnumeric);
            addParameter(p,'fps_variance', 0.7, @isnumeric);
            addParameter(p,'refine_fps', 50, @isnumeric);
            addParameter(p,'C_init', [], @isnumeric);
            parse(p,varargin{:});
            inputs = p.Results;
            settings.max_iter = inputs.max_iter;
            settings.icp_max_iter = inputs.icp_max_iter;
            settings.fps_variance = inputs.fps_variance;
            settings.refine_fps = inputs.refine_fps;
            settings.C_init = inputs.C_init;
        end
        
        function self = set.objectiveFunction_v(self, value)
            assert(isa(value, 'function_handle'), 'objectiveFunction_v has to be a function handle.');
            self.objectiveFunction_v = value;
            self = self.update();
        end
        
        function self = set.objectiveGradient_v(self, value)
            assert(isa(value, 'function_handle'), 'objectiveGradient_v has to be a function handle.');
            self.objectiveGradient_v = value;
            self = self.update();
        end
        
        function self = set.x0_v(self, value)
            assert(isa(value, 'function_handle'), 'x0_v has to be a function handle.');
            self.x0_v = value;
            self = self.update();
        end
        
        function self = set.ref_BASIS_Src(self, value)
            self.ref_BASIS_Src = value;
            self = self.update();
        end
        
        function self = set.ref_BASIS_Tar(self, value)
            self.ref_BASIS_Tar = value;
            self = self.update();
        end
        
        function self = set.est_rank(self, value)
            assert(isnumeric(value), 'est_rank has to be numerical.');
            self.est_rank = value;
            self = self.update();
        end
        
        function self = set.max_iter(self, value)
            assert(isnumeric(value), 'max_iter has to be numerical.');
            self.max_iter = value;
            self = self.update();
        end
        
        function self = set.icp_max_iter(self, value)
            assert(isnumeric(value), 'icp_max_iter has to be numerical.');
            self.icp_max_iter = value;
            self = self.update();
        end
        
        function self = set.fps_variance(self, value)
            assert(isnumeric(value), 'fps_variance has to be numerical.');
            self.fps_variance = value;
            self = self.update();
        end
        
        function self = set.refine_fps(self, value)
            assert(isnumeric(value), 'refine_fps has to be numerical.');
            self.refine_fps = value;
            self = self.update();
        end
    end
    
end

