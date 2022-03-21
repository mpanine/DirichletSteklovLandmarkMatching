function [ P21 ] = tex_transfer_P_Fmap( T12, S1,S2,B1,B2 )
    %TEX_TRANSFER_FMAP 
    B2_inv = pinv(B2);
    B1_inv = pinv(B1);
    % T12: maps S1 -> S2, n1-dim vector
    P21 = sparse(1:S1.nv,double(T12), ones(1,S1.nv), S1.nv, S2.nv); % n1-by-n2 matrix, maps S2 -> S1
    C21 = B1_inv*P21*B2; % k1-by-k2 matrix, maps S2 -> S1
    P21 = B1*C21*B2_inv;
end

