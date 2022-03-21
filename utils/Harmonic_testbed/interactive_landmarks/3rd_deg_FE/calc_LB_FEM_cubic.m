%
% Written by Luca Moschella
% Sapienza University of Rome, 2019
%
function [Stiff, Mass, LumpedMass] = calc_LB_FEM_cubic(M)

Ia = 1/13440 .* ...
    [76 11 11 18 0 27 27 0 18 36;
    11 76 11 0 18 18 0 27 27 36;
    11 11 76 27 27 0 18 18 0 36;
    18 0 27 540 -189 -135 -54 -135 270 162;
    0 18 27 -189 540 270 -135 -54 -135 162;
    27 18 0 -135 270 540 -189 -135 -54 162;
    27 0 18 -54 -135 -189 540 270 -135 162;
    0 27 18 -135 -54 -135 270 540 -189 162;
    18 27 0 270 -135 -54 -135 -189 540 162;
    36 36 36 162 162 162 162 162 162 1944];

Ib = 1/80 .* ...
    [34 -7 0 -54 27 -3 -3 3 3 0;
    -7 34 0 27 -54 3 3 -3 -3 0;
    0 0 0 0 0 0 0 0 0 0;
    -54 27 0 135 -108 0 0 0 0 0;
    27 -54 0 -108 135 0 0 0 0 0;
    -3 3 0 0 0 135 -27 27 27 -162;
    -3 3 0 0 0 -27 135 -135 27 0;
    3 -3 0 0 0 27 -135 135 -27 0;
    3 -3 0 0 0 27 27 -27 135 -162;
    0 0 0 0 0 -162 0 0 -162 324];

Ic = 1/80 .* ...
    [34 0 -7 3 3 -3 -3 27 -54 0;
    0  0 0 0 0 0 0 0 0 0;
    -7 0 34 -3 -3 3  3 -54 27 0;
    3 0 -3 135 -27 27 27 0 0 -162;
    3 0 -3 -27 135 -135 27 0 0 0;
    -3 0 3 27 -135 135 -27 0 0 0;
    -3 0 3 27 27 -27 135 0 0 -162;
    27 0 -54 0 0 0 0 135 -108 0;
    -54 0 27 0 0 0 0 -108 135 0;
    0  0 0 -162 0 0 -162 0 0 324];

Id = 1/80 .* ...
    [68 -7 -7 -51 30 -6 -6 30 -51 0;
    -7 0 7 24 -57 57 -24 0 0 0;
    -7 7 0 0 0 -24 57 -57 24 0;
    -51 24 0 135 -108 27 27 -27 135 -162;
    30 -57 0 -108 135 -135 27 -27 -27 162;
    -6 57 -24 27 -135 135 54 27 27 -162;
    -6 -24 57 27 27 54 135 -135 27 -162;
    30 0 -57 -27 -27 27 -135 135 -108 162;
    -51 0 24 135 -27 27 27 -108 135 -162;
    0 0 0 -162 162 -162 -162 162 -162 324];

q = size(Ib, 1);

Ad = calc_indices_adj_cubic(M);
Nedges = nnz(triu(Ad));
Ntot = M.n + 2*Nedges + M.m;

TRIE = zeros(M.m, 3);
TRIE(:, 1) = Ad(sub2ind(size(Ad), M.TRIV(:, 1), M.TRIV(:, 2))) + M.n + ...
    (M.TRIV(:, 1) > M.TRIV(:, 2));
TRIE(:, 2) = Ad(sub2ind(size(Ad), M.TRIV(:, 1), M.TRIV(:, 2))) + M.n + ...
    (~(M.TRIV(:, 1) > M.TRIV(:, 2)));

TRIE(:, 3) = Ad(sub2ind(size(Ad), M.TRIV(:, 2), M.TRIV(:, 3))) + M.n + ...
    (M.TRIV(:, 2) > M.TRIV(:, 3));
TRIE(:, 4) = Ad(sub2ind(size(Ad), M.TRIV(:, 2), M.TRIV(:, 3))) + M.n + ...
    (~(M.TRIV(:, 2) > M.TRIV(:, 3)));

TRIE(:, 5) = Ad(sub2ind(size(Ad), M.TRIV(:, 3), M.TRIV(:, 1))) + M.n + ...
    (M.TRIV(:, 3) > M.TRIV(:, 1));
TRIE(:, 6) = Ad(sub2ind(size(Ad), M.TRIV(:, 3), M.TRIV(:, 1))) + M.n + ...
    (~(M.TRIV(:, 3) > M.TRIV(:, 1)));

TRIE(:, 7) = (1:M.m)' + M.n + 2*Nedges;
TRItot = [M.TRIV, TRIE];

P1 = M.VERT(M.TRIV(:, 2), :) - M.VERT(M.TRIV(:, 1), :);
P2 = M.VERT(M.TRIV(:, 3), :) - M.VERT(M.TRIV(:, 1), :);

%     Row-wise dot product, repmat to match I#
p11 = dot(P1, P1, 2);
p11b = repmat(p11', q, 1);
p11b = p11b(:);

p22 = dot(P2, P2, 2);
p22b = repmat(p22', q, 1);
p22b = p22b(:);

p12 = dot(P1, P2, 2);
p12b = repmat(p12', q, 1);
p12b = p12b(:);

pre2 = vecnorm(cross(P1,P2, 2), 2, 2);
pre2b = repmat(pre2', q, 1);
pre2b = pre2b(:);

I1b = repmat(Ib, M.m, 1);
I2b = repmat(Ic, M.m, 1);
I3b = repmat(Id, M.m, 1);
I4b = repmat(Ia, M.m, 1);

Alocal2 = (p11b .* I2b + p22b .* I1b - p12b .* I3b) ./ pre2b;
Blocal2 = I4b .* pre2b;

va = Alocal2';
va = va(:);

vb = Blocal2';
vb = vb(:);

idx_rows = repmat(1:M.m, q*q, 1);
idx_rows = idx_rows(:);

idx_cols = repmat(1:q, q, 1);
idx_cols = idx_cols(:);
idx_cols = repmat(idx_cols, M.m, 1);
rows = TRItot(sub2ind(size(TRItot), idx_rows, idx_cols));

idx_cols = repmat(1:q, 1, M.m * q)';
cols = TRItot(sub2ind(size(TRItot), idx_rows, idx_cols));

A = sparse(rows, cols, va, Ntot, Ntot);
B = sparse(rows, cols, vb, Ntot, Ntot);

Stiff = A;
Mass = B;
LumpedMass = sparse(1:M.n, 1:M.n, sum(B, 2), M.n, M.n);
end

function Ad = calc_indices_adj_cubic(M)
    indicesI = [M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,3); M.TRIV(:,2); M.TRIV(:,1)];
    indicesJ = [M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,1); M.TRIV(:,3)];
    Ad = sparse(indicesI, indicesJ, ones(M.m,6), M.n, M.n);
    [indicesI, indicesJ, ~] = find(triu(Ad));
    Nedges = nnz(triu(Ad));
    Ad = sparse(indicesI, indicesJ, 1:2:2*Nedges, M.n, M.n);
    Ad = Ad + Ad';
end
