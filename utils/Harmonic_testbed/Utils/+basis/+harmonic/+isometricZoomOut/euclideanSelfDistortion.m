function dist = euclideanSelfDistortion(S1, p2p)

pre_dist = sum((S1.surface.VERT(p2p,:) - S1.surface.VERT).^2, 2);

dist = sqrt( sum( diag(S1.A).* pre_dist) ) / S1.sqrt_area;



end