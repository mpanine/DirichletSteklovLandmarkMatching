function [w,W] = l1_weight(tri,V,f)
x = V(:,1);
y = V(:,2);
z = V(:,3);
nv = length(x);
nf = size(tri,1);
% generate a function
% f = ones(nv,1);
% f(1) = 1;
% f(2) = -2;
% f(3) = 1;
f_tri = f(tri);
% for each triangle k
%   for each vertex j in triangle k
%     calculate integral over triangle f of 

%calculate triangle areas
L1 = [x(tri(:,2))-x(tri(:,1)),y(tri(:,2))-y(tri(:,1)),z(tri(:,2))-z(tri(:,1))];
L2 = [x(tri(:,3))-x(tri(:,1)),y(tri(:,3))-y(tri(:,1)),z(tri(:,3))-z(tri(:,1))];
Atri = sqrt(sum(cross(L1,L2,2).^2,2))/2;

pos_tri_inds = find(sum(f_tri>=0,2)==3);
neg_tri_inds = find(sum(f_tri<=0,2)==3);
zero_tri_inds = find(sum(f_tri==0,2)==3);
II = kron([pos_tri_inds;neg_tri_inds;zero_tri_inds],[1 1 1]');
JJ = reshape(tri(pos_tri_inds,:)',[],1);
JJ = [JJ;reshape(tri(neg_tri_inds,:)',[],1)];
JJ = [JJ;reshape(tri(zero_tri_inds,:)',[],1)];
SS = kron(Atri(pos_tri_inds)./3,[1 1 1]');
SS = [SS;kron(-Atri(neg_tri_inds)./3,[1 1 1]')];
SS = [SS;kron(zeros(length(zero_tri_inds),1),[1 1 1]')];

W = sparse(II,JJ,SS,nf,nv);
mixed_tri_inds = find(sum(f_tri>0,2)>0 & sum(f_tri<0,2)>0 );

% dW = cell(3,1);
% [dW{:}] = deal(sparse(nf,nv));
% dw_df = sparse(nv,nv);
for k = mixed_tri_inds';
    % find the baricentric coordinates of the zero crossings
    vol_tri = Atri(k)./3;
    curr_tri = tri(k,:);
    f_curr_tri = f_tri(k,:);
%     disp(f_curr_tri);
    zero_ind = find(f_curr_tri==0);
    bc_zero_cross = zeros(2,3);
    if ~isempty(zero_ind)
        bc_zero_cross(1,zero_ind) = 1;
        f_edge = f_curr_tri(1:3~=zero_ind);
        bc_zero_cross(2,(1:3~=zero_ind)) = abs(f_edge)./sum(abs(f_edge));
        
    else
        ind_uniq_sign = find(sign(f_curr_tri)~=sign(sum(sign(f_curr_tri))));
        ind_others = find(1:3~=ind_uniq_sign);
        f_edge = f_curr_tri([ind_uniq_sign,ind_others(1)]);
        bc_zero_cross(1,ind_uniq_sign) = abs(f_edge(2))./sum(abs(f_edge));
        bc_zero_cross(1,ind_others(1)) = abs(f_edge(1))./sum(abs(f_edge));
        f_edge = f_curr_tri([ind_uniq_sign,ind_others(2)]);
        bc_zero_cross(2,ind_uniq_sign) = abs(f_edge(2))./sum(abs(f_edge));
        bc_zero_cross(2,ind_others(2)) = abs(f_edge(1))./sum(abs(f_edge));
    end
    
    
    % calculate the area of the right section
    ind = find(sum(bc_zero_cross>0,1)==2);
    if isempty(ind)
        ind = find(sum(abs(bc_zero_cross-1)<1e-20,1)<1);
        ind = ind(1);
    end;
    
    rem_inds = find(1:3~=ind);
    verts = [x(curr_tri) y(curr_tri) z(curr_tri)];
    v_cross = bc_zero_cross*verts;
    v_spec = verts(ind,:);
    e = verts(rem_inds,:)-repmat(v_spec,2,1);
%     angle_ind = acos(dot(e(1,:),e(2,:))./(norm(e(1,:))*norm(e(2,:))));
    es = v_cross-repmat(v_spec,2,1);
    rem_e = e-es;
%     A_right_section = norm(es(1,:))*norm(es(2,:))*sin(angle_ind)./2;
    A_right_section = norm(cross(es(1,:),es(2,:)))./2;
    % calculate useful derivatives for chain rule later
    
%     des1_dfind = sign(f_curr_tri(ind))*norm(rem_e(1,:))/(abs(f_curr_tri(ind))+abs(f_curr_tri(rem_inds(1))));
%     des2_dfind = sign(f_curr_tri(ind))*norm(rem_e(2,:))/(abs(f_curr_tri(ind))+abs(f_curr_tri(rem_inds(2))));
%     des1_dfrem1 = sign(f_curr_tri(ind))*norm(es(1,:))/(abs(f_curr_tri(ind))+abs(f_curr_tri(rem_inds(1))));
%     des1_dfrem2 = 0;
%     des2_dfrem1 = 0;
%     des2_dfrem2 = sign(f_curr_tri(ind))*norm(es(2,:))/(abs(f_curr_tri(ind))+abs(f_curr_tri(rem_inds(2))));
%     e1 = norm(e(1,:));
%     e2 = norm(e(2,:));
%     f_ind = f_curr_tri(ind);
%     f_ri1 = f_curr_tri(rem_inds(1));
%     f_ri2 = f_curr_tri(rem_inds(2));
    
%     des1_dfind = sign(f_ind)*e1*abs(f_ri1)/((f_ind-f_ri1)^2);
%     des2_dfind = sign(f_ind)*e2*abs(f_ri2)/((f_ind-f_ri2)^2);
%     des1_dfrem1 = e1*f_ind/((f_ind-f_ri1)^2);
%     des1_dfrem2 = 0;
%     des2_dfrem1 = 0;
%     des2_dfrem2 = e2*f_ind/((f_ind-f_ri2)^2);
    
%     dA_dfind = -sign(f_curr_tri(ind))*(norm(es(2,:))*sin(angle_ind)/2*des1_dfind+ norm(es(1,:))*sin(angle_ind)/2*des2_dfind);
%     dA_dfrem1 = -sign(f_curr_tri(ind))*(norm(es(2,:))*sin(angle_ind)/2*des1_dfind);
%     dA_dfrem2 = -sign(f_curr_tri(ind))*(norm(es(1,:))*sin(angle_ind)/2*des2_dfind);
    
%     dA_dfind = (norm(cross(e(1,:),es(2,:)))/2)*des1_dfind/e1 + ...
%                (norm(cross(e(2,:),es(1,:)))/2)*des2_dfind/e2;
%     dA_dfrem1 = (norm(cross(e(1,:),es(2,:)))/2)*des1_dfrem1/e1;
%     dA_dfrem2 = (norm(cross(e(2,:),es(1,:)))/2)*des2_dfrem2/e2;
    
%     dvol_rem_rem_ind1_dfind = (1/(3*norm(e(1,:))))*(des1_dfind*A_right_section + norm(es(1,:))*dA_dfind);
%     dvol_rem_rem_ind2_dfind = (1/(3*norm(e(2,:))))*(des2_dfind*A_right_section + norm(es(2,:))*dA_dfind);
%     dvol_rem_ind_dfind = mean([bc_zero_cross(:,ind);1])*dA_dfind +...
%                          A_right_section/3*(des1_dfind/norm(e(1,:)) + des2_dfind/norm(e(2,:)));
%     dvol_rem_rem_ind1_dfrem1 = (1/(3*norm(e(1,:))))*(des1_dfrem1*A_right_section + norm(es(1,:))*dA_dfrem1);
%     dvol_rem_rem_ind2_dfrem1 = (1/(3*norm(e(2,:))))*norm(es(2,:))*dA_dfrem1;
%     dvol_rem_ind_dfrem1 = mean([bc_zero_cross(:,ind);1])*dA_dfrem1 +...
%                          A_right_section/3*(des1_dfrem1/norm(e(1,:)));
%     dvol_rem_rem_ind1_dfrem2 = (1/(3*norm(e(1,:))))*norm(es(1,:))*dA_dfrem2;
%     dvol_rem_rem_ind2_dfrem2 = (1/(3*norm(e(2,:))))*(des2_dfrem2*A_right_section + norm(es(2,:))*dA_dfrem2);
%     dvol_rem_ind_dfrem2 = mean([bc_zero_cross(:,ind);1])*dA_dfrem1 +...
%                          A_right_section/3*(des1_dfrem1/norm(e(1,:)));
%     dvol_rem_ind_dfind = 2*sign(f_ind)*(mean([bc_zero_cross(:,ind);1])*dA_dfind -...
%                           A_right_section/3*(des1_dfind/e1 + des2_dfind/e2));                  
%     dvol_rem_ind_dfrem1 = 2*sign(f_ind)*(mean([bc_zero_cross(:,ind);1])*dA_dfrem1 -...
%                           A_right_section/3*(des1_dfrem1/e1));                  
%     dvol_rem_ind_dfrem2 = 2*sign(f_ind)*(mean([bc_zero_cross(:,ind);1])*dA_dfrem2 -...
%                           A_right_section/3*(des2_dfrem2/e2));
%     
%     dvol_rem_rem_ind1_dfind = 2*sign(f_ind)*(A_right_section/3*des1_dfind/e1 + dA_dfind*norm(es(1,:))/e1/3);
%     dvol_rem_rem_ind1_dfrem1 = 2*sign(f_ind)*(A_right_section/3*des1_dfrem1/e1 + dA_dfrem1*norm(es(1,:))/e1/3);
%     dvol_rem_rem_ind1_dfrem2 = 2*sign(f_ind)*(dA_dfrem2*norm(es(1,:))/e1/3);
    
    
%     dvol_rem_rem_ind2_dfind = 2*sign(f_ind)*(A_right_section/3*des2_dfind/e2 + dA_dfind*norm(es(2,:))/e2/3);
%     dvol_rem_rem_ind2_dfrem1 = 2*sign(f_ind)*(dA_dfrem1*norm(es(2,:))/e2/3);
%     dvol_rem_rem_ind2_dfrem2 = 2*sign(f_ind)*(A_right_section/3*des2_dfrem2/e2 + dA_dfrem2*norm(es(2,:))/e2/3);
%     
%     dw_df(curr_tri(ind),curr_tri(ind)) = dw_df(curr_tri(ind),curr_tri(ind)) + dvol_rem_ind_dfind;
%     dw_df(curr_tri(rem_inds(1)),curr_tri(ind)) = dw_df(curr_tri(rem_inds(1)),curr_tri(ind)) + dvol_rem_rem_ind1_dfind;
%     dw_df(curr_tri(rem_inds(2)),curr_tri(ind)) = dw_df(curr_tri(rem_inds(2)),curr_tri(ind)) + dvol_rem_rem_ind2_dfind;
%     dw_df(curr_tri(ind),curr_tri(rem_inds(1))) = dw_df(curr_tri(ind),curr_tri(rem_inds(1))) + dvol_rem_ind_dfrem1;
%     dw_df(curr_tri(rem_inds(1)),curr_tri(rem_inds(1))) = dw_df(curr_tri(rem_inds(1)),curr_tri(rem_inds(1))) + dvol_rem_rem_ind1_dfrem1;
%     dw_df(curr_tri(rem_inds(2)),curr_tri(rem_inds(1))) = dw_df(curr_tri(rem_inds(2)),curr_tri(rem_inds(1))) + dvol_rem_rem_ind2_dfrem1;
%     dw_df(curr_tri(ind),curr_tri(rem_inds(2))) = dw_df(curr_tri(ind),curr_tri(rem_inds(2))) + dvol_rem_ind_dfrem2;
%     dw_df(curr_tri(rem_inds(1)),curr_tri(rem_inds(2))) = dw_df(curr_tri(rem_inds(1)),curr_tri(rem_inds(2))) + dvol_rem_rem_ind1_dfrem2;
%     dw_df(curr_tri(rem_inds(2)),curr_tri(rem_inds(2))) = dw_df(curr_tri(rem_inds(2)),curr_tri(rem_inds(2))) + dvol_rem_rem_ind2_dfrem2;
    
    
%     dW{ind}(k,curr_tri(ind)) = 2*sign(f_curr_tri(ind))*dvol_rem_ind_dfind;
%     dW{ind}(k,curr_tri(rem_inds(1))) = 2*sign(f_curr_tri(ind))*dvol_rem_rem_ind1_dfind;
%     dW{ind}(k,curr_tri(rem_inds(2))) = 2*sign(f_curr_tri(ind))*dvol_rem_rem_ind2_dfind;
%     dW{rem_inds(1)}(k,curr_tri(ind)) = 2*sign(f_curr_tri(ind))*dvol_rem_ind_dfrem1;
%     dW{rem_inds(1)}(k,curr_tri(rem_inds(1))) = 2*sign(f_curr_tri(ind))*dvol_rem_rem_ind1_dfrem1;
%     dW{rem_inds(1)}(k,curr_tri(rem_inds(2))) = 2*sign(f_curr_tri(ind))*dvol_rem_rem_ind2_dfrem1;
%     dW{rem_inds(2)}(k,curr_tri(ind)) = 2*sign(f_curr_tri(ind))*dvol_rem_ind_dfrem2;
%     dW{rem_inds(2)}(k,curr_tri(rem_inds(1))) = 2*sign(f_curr_tri(ind))*dvol_rem_rem_ind1_dfrem2;
%     dW{rem_inds(2)}(k,curr_tri(rem_inds(2))) = 2*sign(f_curr_tri(ind))*dvol_rem_rem_ind2_dfrem2;
    % for each vertex in the triangle, calculate theremaining volume
    vol_rem(rem_inds(1)) = A_right_section/3*norm(es(1,:))/norm(e(1,:));
    vol_rem(rem_inds(2)) = A_right_section/3*norm(es(2,:))/norm(e(2,:));
    vol_rem(ind) = sum([bc_zero_cross(:,ind);1])./3*A_right_section;
    
    % find if the base area should be positive or negative
%     if (f_curr_tri(ind)>0) % inside the remaining volume the function is positive
        vol_final = sign(f_curr_tri(ind))*(2*vol_rem - vol_tri);
%         vol_final = (2*vol_rem - vol_tri);
%     else
%         vol_final = vol_tri - 2*vol_rem;
%     end;
    
    % update the weight matrix
    W(k,curr_tri) = vol_final;
end
w = full(sum(W,1));
% W = full(W);
% disp(full(dw_df));
% dw_df = dw_df;
% dw_df = sparse([1:nv,1:nv,1:nv],[tri(:,1)' tri(:,1)' tri(:,1)' tri(:,2)' tri(:,2)' tri(:,2)' tri(:,3)' tri(:,3)' tri(:,3)'],[sum(dW{1}) sum(dW{2}) sum(dW{3})],nv,nv);
