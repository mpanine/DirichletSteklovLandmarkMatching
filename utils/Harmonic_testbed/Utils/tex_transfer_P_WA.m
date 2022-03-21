function [ P21 ] = tex_transfer_P_WA( S1,S2, lmks_S1, lmks_S2 )
    %TEX_TRANSFER_P_WA 
    [E1] = WA_precompute(S1.surface.VERT,S1.surface.TRIV);
    [E2] = WA_precompute(S2.surface.VERT,S2.surface.TRIV);

    A_S1 = vIndex2baryCoords( S1, lmks_S1 );
    A_S2 = vIndex2baryCoords( S2, lmks_S2 );

    pF_12 = compute_pMap( S1, S2, E1, E2, A_S1, A_S2 )';
    
    P21 = sparse(pF_12,1:S1.nv,ones(1,S1.nv),S2.nv,S1.nv);
    
end

