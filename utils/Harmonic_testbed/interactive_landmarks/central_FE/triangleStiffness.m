function stiff = triangleStiffness(triangle, edge_lengths, landmark, edge, radius)
%TRIANGLESTIFFNESS computes the stiffness assosiated with an edge of a SINGLE triangle
%triangle: indices of the vertices
%edge_lengths: lengths of the edges
%landmark: index of the considered landmark (this code does not support having more than one landmark per triangle)
%edge: indices of the edge of the triangle along which the stiffness is computed. Unlike actual edges, can be the same index repeated twice.
%radius: radius of the central disk at the landmark


[l1, tmax, a, b, Rloc, ind1, ind2] = triangleStandardForm(triangle, edge_lengths, landmark, edge, radius ); %Standard form



% [f1t,f1r] = gradU(ind1, l1, a, b, Rloc); %Functions to integrate    %%OLD VERSION
% [f2t,f2r] = gradU(ind2, l1, a, b, Rloc);
% 
% 
% 
% stiff = centralElementIntegral(f1t, f1r, f2t, f2r, tmax, Rloc, maxR(l1, a, b));


f = stiffnessIntegrand(l1, a, b, Rloc, ind1, ind2);

stiff = integral(f, 0 ,tmax);

% th = linspace(0,tmax,1000);
% plot(th,f(th))
% hold on

end