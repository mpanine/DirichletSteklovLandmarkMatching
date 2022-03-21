function [obj_fun,const_viol] = compressed_modes_objective(tri,V,W,A,mu,F)

obj_fun = trace(F'*W*F) + mu*l1_fnorm(tri,V,F);
disp(['Obj func: ' num2str(obj_fun)]);
const_viol = norm(F'*A*F - speye(size(F,2)));
disp(['Const violation: ' num2str(const_viol)]);