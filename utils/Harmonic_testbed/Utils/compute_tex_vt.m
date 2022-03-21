function [ S_vt ] = compute_tex_vt( S,scale )
    %COMPUTE_TEX_VT 
    cols = [1 2 3];
    X = S.surface.VERT;
    [~,ic] = sort(std(X),'descend');
    S_vt = shape.io.generate_tex_coords(X(:,ic), cols(1), cols(2), scale);
    
end

