function [vs,ds]=hamiltonian_basis(Q,B,card,iso)

%ver: matrix of vertex coordinates n*3
%tri: matrix of triangles t*3
%index: index of the vertex for the neighbor computation
%depth: depth of the neighbor of the vertex 'index'
%Q: FEM stiffness matrix
%B: FEM mass matrix
%card: number of the computed Hamiltonian eigenvalues and eigenvectors
%iso: value of the potential (e.g., iso=100.0);
%ds: matrix of Hamiltonian eigenvalues
%vs: matrix of Hamiltonian eigenvectors

% num_ver=size(ver,1); num_tri=size(tri,1);
% next=[2 3 1]; vv=sparse(num_ver,num_ver);
% for i=1:3
%     vv=vv+sparse(tri(:,i),tri(:,next(i)),1,num_ver,num_ver,num_tri);
% end
% 
% aux=[index];
% for i=1:depth
%     new_aux=[];
%     for j=1:length(aux)
%         tmp=find( vv(aux(j),:) );
%         new_aux=[new_aux,tmp];
%     end
%     aux=unique(new_aux);
% end
% %aux list of indices of the local neighbor at vertex index
% list=[1:num_ver]; list(aux)=[];
% 
% V=sparse(list,list,1.0,num_ver,num_ver,length(list));

L=-Q+B*diag(iso);

[vs,ds]=eigs(L,B,card,'sm');
[ds,sort_idx] = sort(diag(ds));
vs = vs(:,sort_idx);