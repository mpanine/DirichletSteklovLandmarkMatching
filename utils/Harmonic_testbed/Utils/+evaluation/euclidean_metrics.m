function [ metrics ] = euclidean_metrics( Tar, pF, indices_Src, indices_Tar )
%EUCLIDEAN_METRICS Computes the euclidean error made by the input p2p map at the Source and Target indices.
    val = sqrt(sum(...
    (Tar.surface.VERT(indices_Tar,:) - ...
     Tar.surface.VERT(pF(indices_Src),:)).^2,2))/Tar.sqrt_area;

    metrics.values = val;
    metrics.mean = mean(val);
    metrics.worstCase = max(max(val));
    metrics.standardDeviation = std(val);

end

