function [ metrics ] = geodesic_metrics( Tar, pF, indices_Src, indices_Tar )
%GEODESIC_METRICS Computes the geodesic error made by the input p2p map at the Source and Target indices.
    assert(length(indices_Src)==length(indices_Tar),...
            'GEODESICMETRICSERROR:lengthUnequal', 'The number of indices on Source and Target shapes have to be equal.');
        
    pF_indices_Src = pF(indices_Src);
    if(size(indices_Tar, 2) > size(indices_Tar,1))
            indices_Tar = indices_Tar';
    end
    if(size(pF_indices_Src, 2) > size(pF_indices_Src,1))
            pF_indices_Src = pF_indices_Src';
    end
        
    if ~isfield(Tar,'Gamma')
        val = shape.geometry.geodesics_pairs(Tar,[indices_Tar,pF_indices_Src])/Tar.sqrt_area;
    else
        samples = 1:Tar.nv;
        Gamma = Tar.Gamma;
        
        max_geo_dist = max(max(Gamma(~isinf(Gamma))));
        Gamma(isinf(Gamma)) = max_geo_dist;

        [~, samples_Tar] = ismember(indices_Tar, samples);
        [~, samples_Src] = ismember(pF_indices_Src, 1:Tar.nv);
        val = diag(Gamma(samples_Tar,samples_Src)/Tar.sqrt_area);
    end
    metrics.values = val;
    metrics.mean = mean(val);
    metrics.worstCase = max(max(val));
    metrics.standardDeviation = std(val);

end

