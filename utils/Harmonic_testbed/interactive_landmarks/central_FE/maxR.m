function mR = maxR(l1, a, b)
%MAXR produces a function that gives the maximum radius (opposite side of
%triangle) as a function of 



mR = @(t) b.*l1.*(b.*cos(t)+((-1).*a+l1).*sin(t)).^(-1);



end

