function hks = heatKernelSignature(laplaceBasis, eigenvalues, Ae, numTimes)
% This method computes the heat kernel signature for each vertex on a list.
% It uses precomputed LB eigenstuff stored in "mesh" and automatically
% chooses the time steps based on mesh geometry.

numVertices = size(laplaceBasis, 1);

tmin = -4*log(10)/eigenvalues(end);
tmax = -4*log(10)/eigenvalues(10); % Why 10?
ts = logspace(log10(tmin),log10(tmax),numTimes);

areaMtx = sparse(Ae);
D = laplaceBasis' * (areaMtx * laplaceBasis.^2);

T = exp(-abs(eigenvalues)*ts);
F = D*T;

fs = laplaceBasis * F;    
nf = size(F,2);
scale = sparse(1:nf, 1:nf, 1./max(fs), nf, nf);
fs = fs*scale;    
hks =  laplaceBasis'*(areaMtx*fs);

hks = real(laplaceBasis * hks);