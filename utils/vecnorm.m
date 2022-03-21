function n = vecnorm(A, p, dim)
if(nargin < 2)
    p = 2; dim = find(size(A)>1,1);
elseif (nargin < 3) 
    dim = find(size(A)>1,1);
end
A2p = A.^p;
s   = sum(A2p, dim);
n   = s.^(1/p);