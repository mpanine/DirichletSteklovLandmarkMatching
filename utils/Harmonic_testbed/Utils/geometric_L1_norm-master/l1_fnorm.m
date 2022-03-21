function l1_nrm = l1_fnorm(tri,V,F)
l1_nrm = zeros(1,size(F,2));
for ii = 1:size(F,2);
    [w] = l1_weight(tri,V,F(:,ii));
    l1_nrm(ii) = sum(w'.*F(:,ii));
end;