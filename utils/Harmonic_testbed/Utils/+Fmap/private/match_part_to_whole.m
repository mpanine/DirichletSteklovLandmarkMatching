%
% This code accompanies the paper:
%
% "Partial Functional Correspondence"
% Rodola, Cosmo, Bronstein, Torsello, Cremers
% Computer Graphics Forum 2016
%
% Please cite the paper above if you use this code in your research.
%
% Written by Emanuele Rodola and Luca Cosmo
%
function [C, v, matches] = match_part_to_whole(Src, Tar, G, F, C_init, options)
mu1 = options.mu1;
mu2 = options.mu2;
mu3 = options.mu3;
mu4 = options.mu4;

v_init = ones(Src.SHAPE.nv,1);

C_prev = [];
v_prev = [];
matches_prev = [];
objs = [];

k = size(C_init,1);

for i=1:options.max_iters
    
    fprintf('------------------------- Iteration %d -------------------------\n',i);
    
    [C, est_rank, cost] = optimize_F(Src, Tar, G, F, v_init, C_init, mu1, mu2);
    
    if i>1 && cost >= objs(end)
        fprintf('Cost increasing, stopping now and keeping previous solution.\n');
        C = C_prev;
        v = v_prev;
        matches = matches_prev;
        break;
    end
        
    objs = [objs cost];
    
    [Co, matches] = icp_refine_partial(Src.BASIS, Tar.BASIS, est_rank, C, options.icp_max_iters);
    
    C_init = [Co zeros(size(C,1),size(C,1)-size(Co,2))];

    v = optimize_v(Src, Tar, G, F, C_init, mu3, mu4, options);
    v = (0.5*tanh(6*(v-0.5))+0.5);
    v_init = v;
    
    C_prev = C_init;
    v_prev = v_init;
    matches_prev = matches;
    
end

if i==options.max_iters
    fprintf('Maximum number of iterations reached in the alternating process.\n');
end

fprintf('Done.\n');

end
