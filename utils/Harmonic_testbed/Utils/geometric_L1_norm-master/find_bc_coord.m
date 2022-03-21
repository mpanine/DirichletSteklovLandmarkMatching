function [bc_coord,err] = find_bc_coord(P1,P2,P3,P,flag)
if nargin<5 
  flag = 1;
end;
A =  ([[P1' P2' P3']; 1e3*ones(1,3)]);
b = ([P';1e3]);

bc_coord =A\b;
err = norm(A*bc_coord - b);
if (flag)
if abs(sum(abs(bc_coord))-1)>1e-4 | norm([P1' P2' P3']*bc_coord - P')>1e-4
  bc_coord = nan;
end;
end