function Fmap = CentralSteklovFmapFromP2P(Src_evecs_bound, Tar_evecs_bound, Tar_Mass_Boundary, p2p, ZOsize)
%CentralSteklovFmapFromP2P: extracts functional maps from p2p maps for the central matching.

%p2p: Tar -> Src

Fmap = zeros(ZOsize, ZOsize);

for i = 1:length(Src_evecs_bound)
    

    Fmap = Fmap + Tar_evecs_bound{i}(:,1:ZOsize)' * Tar_Mass_Boundary{i} * Src_evecs_bound{i}(p2p{i},1:ZOsize);


end



end