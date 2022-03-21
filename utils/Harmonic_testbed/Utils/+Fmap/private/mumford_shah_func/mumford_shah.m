%
% This code accompanies the paper:
%
% "Partial Functional Correspondence"
% Rodola, Cosmo, Bronstein, Torsello, Cremers
% Computer Graphics Forum 2016
%
% Please cite the paper above if you use this code in your research.
%
% Written by Luca Cosmo
%
classdef mumford_shah < handle
    properties (Hidden = true, SetAccess = private)
        cpp_handle;
        M_triangles;
        M_areas_onering;
    end
    methods
	
        % Constructor
        function this = mumford_shah(VERT, TRIV, areas_onering)
            this.cpp_handle = mumford_shah_wrapper(0, VERT', TRIV');
            this.M_triangles = TRIV;
            this.M_areas_onering = areas_onering;
        end
		
        % Destructor
        function delete(this)
            mumford_shah_wrapper(-1, this.cpp_handle);
        end

        function val = cost(this, CA, PtS, G, v, perturb, target_area, areas, mu1, mu2, options)
            tanhv = tanh(6*(v-0.5));
            tv = (0.5*tanhv+0.5);
            VG = bsxfun(@times,G,tv);
            
            data_term = sum( sqrt(sum((CA - PtS * VG).^2,1)));
			area_term = (target_area - sum((0.5*tanhv+0.5).*areas)).^2;
            
            reg_term = sum( mumford_shah_wrapper(1,this.cpp_handle, v, options.tv_sigma, 0, options.tv_mean, perturb) );
            
            val = data_term + mu1*area_term + mu2*reg_term;		
        end
        
        function val = grad(this, CA, PtS, G, v, perturb, target_area, areas, mu1, mu2, options)
            tanhv = tanh(6*(v-0.5));
            B = full(PtS * diag(sparse((0.5*tanhv+0.5))) * G);

            data_term = gradient_v_L21(CA-B,1./sqrt(sum((CA-B).^2)),PtS,3*(1-tanhv.^2),G);
          
			area_term = -6*(target_area - sum((0.5*tanhv+0.5).*areas))*(1 - tanhv.^2).*areas;

            reg_term = mumford_shah_wrapper(2,this.cpp_handle, v, options.tv_sigma, 0, options.tv_mean, perturb);

            val = data_term + mu1*area_term + mu2*reg_term;		
        end
    end
end