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
function [Fs] = compute_indicator_functions(Ms,matches,alpha)

	% NOTE: fastmarchmex computes approximate geodesics using fast marching.
	% You can use your own code for the same task here.
	
    ncorr = size(matches,2);
	
    parfor m = 1:numel(Ms)
	
        M = Ms{m};
        n = size(M.surface.VERT,1); 
        F = [];
		
        Mdist = fastmarchmex('init', int32(M.surface.TRIV-1), double(M.surface.VERT(:,1)), double(M.surface.VERT(:,2)), double(M.surface.VERT(:,3)));
		
        for i=1:ncorr
            if(matches(m,i)<=0)
                F(:,i) = zeros(n,1);
                continue;
            end
            source = inf(n, 1);
            source(matches(m,i)) = 0;
            dM = exp(-0.5*alpha*fastmarchmex('march', Mdist, double(source)));
            F(:,i)=dM;   
        end
		
        Fs{m}=F;
        fastmarchmex('deinit', Mdist);
		
    end
	
end
