% Euclidean farthest point sampling.
function idx = euclidean_fps(S,k,seed)

C = S.surface.VERT;
nv = size(C,1);

if(nargin<3)
    idx = randi(nv,1);
else
    idx = seed;
end

dists = bsxfun(@minus,C,C(idx,:));
dists = sum(dists.^2,2);

for i = 1:k
    maxi = find(dists == max(dists));
    maxi = maxi(1);
    idx = [idx; maxi];
    newdists = bsxfun(@minus,C,C(maxi,:));
    newdists = sum(newdists.^2,2);
    dists = min(dists,newdists);    
end

idx = idx(2:end);