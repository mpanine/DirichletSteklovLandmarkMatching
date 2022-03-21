% Compute the Steklov boundary operator (a mass matrix concentrated on the boundary)
function S = boundarySteklov(X, E)
nv = size(X,1);
ne = size(E,1);



%% Non-lumped version -- Sometimes fails to normalize

% L6 = normv(X(E(:,1),:)-X(E(:,2),:))/6; % Sixth of the edge lengths
% 
% S = sparse(E(:,1),E(:,2), L6, nv, nv);
% S = (S + S')/2;
% S = S - sparse(diag(diag(S)));
% 
% d = sum(S);
% 
% S = S + sparse(diag(d));

%% Lumped version 

L = normv(X(E(:,1),:)-X(E(:,2),:)); % Edge lengths
LL = [L; L]; %Two copies of edge lengths, one for each edge endpoint.

relevant_verts = [E(:,1); E(:,2)]; %All the edge endpoints in the same order.

S = sparse(relevant_verts, relevant_verts, LL, nv, nv) / 2; %Half the sum of adjacent edge lengths.



end