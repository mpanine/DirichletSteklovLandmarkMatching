function Fmap = SteklovDirichletFmapFromP2P(Src_evecs_bound, Tar_evecs_bound, Tar_Mass_Boundary, p2p, ZOsize)
%SteklovDirichletFmapFromP2P extracts functional maps from p2p maps.

%p2p: Tar -> Src

Fmap = cell(length(Src_evecs_bound),1);

for i = 1:length(Src_evecs_bound)
    

    Fmap{i} = Tar_evecs_bound{i}(:,1:ZOsize)' * Tar_Mass_Boundary{i} * Src_evecs_bound{i}(p2p{i},1:ZOsize);


end


end