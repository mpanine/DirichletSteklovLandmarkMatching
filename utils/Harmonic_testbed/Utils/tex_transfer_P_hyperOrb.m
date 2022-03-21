function [ P21 ] = tex_transfer_P_hyperOrb( S1,S2, lmks_S1, lmks_S2 )
    %TEX_TRANSFER_P_HYPERORB Summary of this function goes here
    %   Detailed explanation goes here
    embedder_S1=Flattener(S1.surface.VERT,S1.surface.TRIV,lmks_S1);
    embedder_S2=Flattener(S2.surface.VERT,S2.surface.TRIV,lmks_S2);

    flatteners={embedder_S1,embedder_S2};
    mapper=Mapper(flatteners);

    mapper.computeMap();
    
    mapper.lift(2,1);
    P21 = mapper.map.barCoords{1,2};
    
end

