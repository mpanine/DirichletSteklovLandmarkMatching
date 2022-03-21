function [Afun] = Afunc(tri,V,W,A,epsilon,mu,f)
ww = l1_weight(tri,V,f);
ww = sparse(ww);
ww = (abs(ww')./(sqrt(epsilon+f.^2).*A));
ww = diag(ww);
Afun = ((-W+mu*ww)*diag(1./(A)))\f;
end
