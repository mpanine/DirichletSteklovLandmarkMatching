function [ pF_refined_BCICP, pF_inverse_refined_BCICP ] = bcicp( Src, Tar, pF, pF_inverse, numIter )
%COMPUTE_BCICP computes the BCICP refinement
    assert(isfield(Src, 'SHAPE'), 'Source struct requires a SHAPE field containing its shape.');
    assert(isfield(Src, 'BASIS'), 'Source struct requires a BASIS field containing its basis.');
    assert(isfield(Tar, 'SHAPE'), 'Target struct requires a SHAPE field containing its shape.');
    assert(isfield(Tar, 'BASIS'), 'Target struct requires a BASIS field containing its basis.');
    [ pF_refined_BCICP, pF_inverse_refined_BCICP ] = bcicp_refine(Src.SHAPE, Tar.SHAPE, Src.BASIS.basis, Tar.BASIS.basis, pF, pF_inverse, numIter);
end

