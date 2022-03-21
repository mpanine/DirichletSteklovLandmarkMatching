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
classdef mumford_shah_Cu < handle
    properties (Hidden = true, SetAccess = private)
        shapeM;
        shapeN;
        handle_F;
        handle_P;
        calln;
        target_area;
    end
    methods
	
        % Constructor
        function this = mumford_shah_Cu(F,P)
            this.handle_F = mumford_shah_wrapper2(0, F.VERT', F.TRIV');
            this.handle_P = mumford_shah_wrapper2(0, P.VERT', P.TRIV');
            this.calln=0;
            this.shapeM = F;
            this.shapeN = P;
            
            this.target_area= full(sum(diag(P.S)));
        end
		
        % Destructor
        function delete(this)
            mumford_shah_wrapper2(-1, this.handle_F);
            mumford_shah_wrapper2(-1, this.handle_P);
        end

        function val = cost(this, C, Psi, PhitS,options)
            
            
            data_term = 0;
            v = Psi*C*PhitS;
            reg_term1 = (mumford_shah_wrapper2(1,this.handle_F, v, options.tv_sigma, 0, 0, 0));
            reg_term2 = (mumford_shah_wrapper2(1,this.handle_F, v, options.tv_sigma, 0, 1, 0));

            val = data_term - sum(reg_term1) - sum(reg_term2);

        end
        
        function val = grad(this, C, PhitS,options)
            data_term =  0;
            v = this.shapeM.evecs*C*PhitS;
            reg_term1 = -mumford_shah_wrapper2(2,this.handle_F, v, PhitS,this.shapeM.evecs, options.tv_sigma, 0);
            reg_term2 = -mumford_shah_wrapper2(2,this.handle_F, v, PhitS,this.shapeM.evecs, options.tv_sigma, 1);
            val = data_term + reg_term1*PhitS'+ reg_term2*PhitS';	
        end
    end
end