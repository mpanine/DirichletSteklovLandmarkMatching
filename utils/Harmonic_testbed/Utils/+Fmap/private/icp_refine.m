function [C_out, matches] = icp_refine(BASIS_Src, BASIS_Tar, C, nk)
    Vs = BASIS_Src;
    Vt = BASIS_Tar;
    n1 = size(BASIS_Src.basis,2);
    n2 = size(BASIS_Tar.basis,2);

    C_out = C;
    
    for k=1:nk
        matches = knnsearch((C_out*Vs.basis')',Vt.basis);
        W = BASIS_Tar.basis\BASIS_Src.basis(matches,:);
        
        [s,~,d] = svd(W);
        C_out = s*eye(n2,n1)*d';
    end
end