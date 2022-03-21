function [l1_nrm,A] = l1_anorm(tri,V,F)

G=metric_scale(V(:,1),V(:,2),V(:,3),tri,0);
[mat_struct]=compute_div_grad_matrices(V(:,1),V(:,2),V(:,3),tri,G);
A = mat_struct.A./3;
l1_nrm = zeros(1,size(F,2));
for ii = 1:size(F,2);
  l1_nrm(ii) = sum(abs(A.*F(:,ii)));
end;