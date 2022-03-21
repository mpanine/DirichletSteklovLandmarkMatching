classdef Settings < core.Settings
    %SETTINGS class of shape settings
    properties
        computeNormals % compute the face and vertex normals (normals_vtx, normals_face)
        computeMHB % compute the Laplace-Beltrami eigenfunction basis
        MHB_settings % settings of the Laplace-Beltrami eigenfunction basis
        computeGeoDist % find the geodesic distance matrix (Gamma)
        computeObj % compute the .obj file associated to the shape
        findEdge % find the edge list and the corresponding cotangent weight (Elist, EdgeWeight)
        findNeigh % find the one-ring neighbor for each vertex (vtx_neigh)
    end
    
    methods
        function settings = Settings(varargin)
            p = inputParser;
            addParameter(p,'computeNormals',false,@islogical);
            addParameter(p,'computeMHB',false,@islogical);
            addParameter(p,'MHB_settings',basis.MHB_Settings(200));
            addParameter(p,'computeGeoDist',false,@islogical);
            addParameter(p,'computeObj',false,@islogical);
            addParameter(p,'findEdge',true,@islogical);
            addParameter(p,'findNeigh',true,@islogical);
            parse(p,varargin{:});
            inputs = p.Results;
            
            assert(isa(inputs.MHB_settings, 'basis.MHB_Settings'),...
                      sprintf('MHB_settings has to be a basis.MHB_Settings object got %s object instead.',class(inputs.MHB_settings)));
            
            settings.computeNormals = inputs.computeNormals;
            settings.computeMHB = inputs.computeMHB;
            settings.MHB_settings = inputs.MHB_settings;
            settings.computeGeoDist = inputs.computeGeoDist;
            settings.computeObj = inputs.computeObj;
            settings.findEdge = inputs.findEdge;
            settings.findNeigh = inputs.findNeigh;
        end
        
        function self = set.computeNormals(self, value)
            assert(islogical(value), 'computeNormals has to be a logical.');
            self.computeNormals = value;
            self = self.update();
        end
        
        function self = set.computeMHB(self, value)
            assert(islogical(value), 'computeMHB has to be a logical.');
            self.computeMHB = value;
            self = self.update();
        end
        
        function self = set.computeGeoDist(self, value)
            assert(islogical(value), 'computeGeoDist has to be a logical.');
            self.computeGeoDist = value;
            self = self.update();
        end
        
        function self = set.computeObj(self, value)
            assert(islogical(value), 'computeObj has to be a logical.');
            self.computeObj = value;
            self = self.update();
        end
        
        function self = set.findEdge(self, value)
            assert(islogical(value), 'findEdge has to be a logical.');
            self.findEdge = value;
            self = self.update();
        end
        
        function self = set.findNeigh(self, value)
            assert(islogical(value), 'findNeigh has to be a logical.');
            self.findNeigh = value;
            self = self.update();
        end
        
        function self = set.MHB_settings(self, value)
            assert(isa(value, 'basis.MHB_Settings'),...
                      sprintf('MHB_settings has to be a basis.MHB_Settings object got %s object instead.',class(value)));
            self.MHB_settings = value;
            self = self.update();
        end
    end
    
end

