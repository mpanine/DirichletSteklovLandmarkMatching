function [ww] = l1_weight_vec(tri,V,F)
ww = zeros(size(F'));
for ii = 1:size(F,2)
    f = F(:,ii);
    x = V(:,1);
    y = V(:,2);
    z = V(:,3);
    nv = length(x);
    nf = size(tri,1);
    f_tri = f(tri);
    
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
    
    for k = mixed_tri_inds';
        vol_tri = Atri(k)./3;
        curr_tri = tri(k,:);
        f_curr_tri = f_tri(k,:);
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
        es = v_cross-repmat(v_spec,2,1);
        rem_e = e-es;
        A_right_section = norm(cross(es(1,:),es(2,:)))./2;
        vol_rem(rem_inds(1)) = A_right_section/3*norm(es(1,:))/norm(e(1,:));
        vol_rem(rem_inds(2)) = A_right_section/3*norm(es(2,:))/norm(e(2,:));
        vol_rem(ind) = sum([bc_zero_cross(:,ind);1])./3*A_right_section;
        
        vol_final = sign(f_curr_tri(ind))*(2*vol_rem - vol_tri);
        
        W(k,curr_tri) = vol_final;
    end
    w = full(sum(W,1));
    ww(:,ii) = w;
    
end;