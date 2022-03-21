%
% Apply rank-k Sherman-Morrison update to solve (A+eta*B*B')*x=b,
% where k=size(B,2). Very efficient if k<<size(B,1).
%
% Written by Emanuele Rodola
%
function x = sherman(b, A, B, eta)

k = size(B,2);

y = A\b;
z = A\(eta*B);
x = y - z * ( (speye(k)+B'*z) \ (B'*y) );

end
