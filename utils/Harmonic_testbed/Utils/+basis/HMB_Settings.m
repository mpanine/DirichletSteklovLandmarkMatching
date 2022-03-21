classdef HMB_Settings < basis.Settings
    %MHB_SETTINGS Settings for the MHB basis
    
    properties
    end
    
    methods
        function settings = HMB_Settings(varargin)
            p = inputParser;
            addOptional(p,'numBasisFun',1);
            addOptional(p,'numControlLandmarks',0);
            addOptional(p,'intrinsic',false);
            addOptional(p,'numLBFun',50);
            parse(p,varargin{:});
            
            numBasisFun = p.Results.numBasisFun;
            parameters = struct;
            parameters.numControlLandmarks = p.Results.numControlLandmarks;
            parameters.intrinsic = p.Results.intrinsic;
            parameters.numLBFun = p.Results.numLBFun;
            
            settings@basis.Settings(numBasisFun, parameters);
            
            settings.name = 'HMB';
            
            settings = settings.update();
        end
    end
    
end