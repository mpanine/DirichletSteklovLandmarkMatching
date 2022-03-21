function [l1_nrm] = l1_norm(tri,V,F)

G=metric_scale(V(:,1),V(:,2),V(:,3),tri,0);
[mat_struct]=compute_div_grad_matrices(V(:,1),V(:,2),V(:,3),tri,G);
A = mat_struct.A;
n = size(V,1);
tot_A = sum(sum(A));
l1_nrm = zeros(1,size(F,2));
for ii = 1:size(F,2);
  l1_nrm(ii) = tot_A*sum(abs(F(:,ii)))./n;
end;