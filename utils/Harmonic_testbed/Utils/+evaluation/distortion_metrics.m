function [ metrics ] = distortion_metrics( Src, Tar, pF )
%DISTORTION_METRICS Computes the distortion made by the input p2p map at all locations on the target shape.
    assert(isfield(Tar, 'Gamma'),...
    'The Source shape needs to have a Gamma field to compute the distortion metrics. Set computeGeoDist to true when creating the shape.');
    assert(isfield(Tar, 'Gamma'),...
    'The Target shape needs to have a Gamma field to compute the distortion metrics. Set computeGeoDist to true when creating the shape.');
    
    D = Tar.Gamma/max(max(Tar.Gamma))-Src.Gamma(pF,pF)/max(max(Src.Gamma));
    A_x_Tar = repmat(sum(Tar.A,1),[size(Tar.A,1),1]);
    A_y_Tar = repmat(sum(Tar.A,2),[1,size(Tar.A,2)]);
    
    metrics.mean = full((1/sum(sum(Tar.A)))*sqrt(sum(sum((D.^2).*A_x_Tar.*A_y_Tar))));
%     metrics.values = D;
%     metrics.values = [];
    metrics.values = mean(D,1);
    metrics.standardDeviation = full(sqrt(1/(sum(sum(Tar.A)))*sqrt(sum(sum((((D.^2)-metrics.mean*ones(size(D))).^2).*A_x_Tar.*A_y_Tar)))));
    metrics.worstCase = full(max(max(abs(D))));

end

