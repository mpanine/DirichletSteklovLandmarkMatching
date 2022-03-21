clear all;
close all;

addpath(genpath('./utils'));

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings();
%% paths
data_dir = './data/FAUST/';
shape_dir = [data_dir,'vtx_5k/'];
lmk_dir = [data_dir,'Landmarks/'];
geo_dir = [data_dir,'Gammas/'];
if ~exist(lmk_dir,'dir')
    mkdir(lmk_dir);
end
if ~exist(geo_dir,'dir')
    mkdir(geo_dir);
end

save_dir = './data/FAUST/remeshed_mat/';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% Generate shapes
all_shapes = dir([shape_dir,'*.off']);
all_shapes = {all_shapes.name};

r_ls = [0.25,0.5,0.75];

for shape_id=1:size(all_shapes,2)
    shape_name = all_shapes{shape_id};
    shape_name = strrep(shape_name,'.off','');
    if strcmp(shape_name,'tr_reg_002')||strcmp(shape_name,'tr_reg_007')
        SHAPE = shape.compute([shape_dir,shape_name,'.off'],shape_settings,cache_settings);
        load([lmk_dir,shape_name,'.mat']);
        if ~exist([geo_dir,shape_name,'.mat'],'file')
            Gamma = compute_geodesic_dist_matrix(SHAPE,1:size(SHAPE.surface.VERT,1));
            save([geo_dir,shape_name,'.mat'],'Gamma');
        else
            load([geo_dir,shape_name,'.mat']);
        end
        save_struct = struct;
        save_struct.SHAPE = SHAPE;
        save_struct.landmarks = landmarks;
        save_struct.Gamma = Gamma;
        save([save_dir,shape_name,sprintf('_R_%.2f.mat',1.)],'-struct','save_struct');

        for r=r_ls
            [TRI_remeshed,VERT_remeshed] = reducepatch(SHAPE.surface.TRIV,SHAPE.surface.VERT,r);
            shape_re = struct;%a mesh structure with surface (X,Y,Z, TRIV)
            shape_re.surface.X = VERT_remeshed(:,1);
            shape_re.surface.Y = VERT_remeshed(:,2);
            shape_re.surface.Z = VERT_remeshed(:,3);
            shape_re.surface.VERT = VERT_remeshed;
            shape_re.surface.TRIV = TRI_remeshed;
            shape_re.name = [shape_name,sprintf('_re%.2f',r)];
            SHAPE_re = shape.compute(shape_re,shape_settings,cache_settings);

            landmarks_re = knnsearch(SHAPE_re.surface.VERT,SHAPE.surface.VERT(landmarks,:));
            reverse_corres = knnsearch(SHAPE.surface.VERT,SHAPE_re.surface.VERT); % to get the vertex id of the full shape corresponding to the remeshed shape (used to extract reduced Gamma)
            Gamma_re = Gamma(reverse_corres,reverse_corres);

            save_struct = struct;
            save_struct.SHAPE = SHAPE_re;
            save_struct.landmarks = landmarks_re;
            save_struct.Gamma = Gamma_re;
            save_struct.reverse_corres = reverse_corres;
            save([save_dir,shape_name,sprintf('_R_%.2f.mat',r)],'-struct','save_struct');
        end
    end
end