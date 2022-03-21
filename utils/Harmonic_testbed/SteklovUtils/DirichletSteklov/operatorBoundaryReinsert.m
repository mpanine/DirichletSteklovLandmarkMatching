function OP_out = operatorBoundaryReinsert(OP, landmarks)
%OP is a square operator that does not include landmark vertices.
%Reinserts zeros at the landmarks of OP.

ss = length(landmarks) + size(OP,1);

% OP_out = zeros(ss,ss);
OP_out = sparse(ss,ss);


interior = 1:ss;
interior(landmarks) = [];


OP_out(interior,interior) = OP;







end