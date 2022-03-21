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
function [C_ref, v_ref, matches_ref] = refine_partial(Src, Tar, C, matches, options)

% Refine by using landmarks

fps = shape.geometry.euclidean_fps(Tar.SHAPE, options.refine_fps, 1);
FG = compute_indicator_functions({Tar.SHAPE,Src.SHAPE}, [fps matches(fps)]', options.fps_variance);
F = FG{1};
G = FG{2};

[C_ref, v_ref, matches_ref] = ...
    match_part_to_whole(Src, Tar, G, F, C, options);

end
