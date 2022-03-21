function [  ] = tex_transfer_save_models( S1,S2,S1_vt,S2_vt, method,save_path )
    %TEX_TRANSFER_SAVE_MODELS 
    
    texture_im = 'texture.jpg';

    if nargin < 5, save_path = './'; end
    
    S1_name = 'target';
    S2_name = 'source';
    
    if isfield(S1,'name')
        S1_name = [S1_name,'_', S1.name];
    end

    if isfield(S2,'name')
        S2_name = [S2_name,'_', S2.name];
    end

    % Source mesh
    shape.io.wobj_with_texture(S2, S2_vt, texture_im, [save_path, S2_name,'_',method]);
    % Target mesh: 
    shape.io.wobj_with_texture(S1, S1_vt, texture_im, [save_path,S1_name,'_',method]);
    
end

