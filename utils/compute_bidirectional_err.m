function [err1,curve1,err2,curve2] = compute_bidirectional_err(Src,Tar,eval_Src,eval_Tar,p2pSrcTar,p2pTarSrc,thresh)
%COMPUTE_BIDIRECTIONAL_ERR Summary of this function goes here
%   Detailed explanation goes here

    err1 = calc_geo_err(p2pSrcTar(eval_Src),eval_Tar,Tar.gamma);
    curve1 = calc_err_curve(err1,thresh);
    
    err2 = calc_geo_err(p2pTarSrc(eval_Tar),eval_Src,Src.gamma);
    curve2 = calc_err_curve(err2,thresh);
    
    
end

