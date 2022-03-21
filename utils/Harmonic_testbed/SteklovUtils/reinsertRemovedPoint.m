function funcs_re = reinsertRemovedPoint(funcs, point)
% Re-inserts a 0 value into funcs at a SINGLE REMOVED POINT

ss = size(funcs);

funcs_re = zeros(ss(1) + 1 , ss(2));

if point > 1
    funcs_re(1:point - 1,:) = funcs(1:point - 1,:);
end

if point < size(funcs,1)
    funcs_re(point + 1:end,:) = funcs(point:end,:);
end



end