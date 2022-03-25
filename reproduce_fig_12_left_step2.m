clear all;
close all;

addpath(genpath('./utils'));

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings();

data_dir = './data/FAUST/';
mat_dir = [data_dir,'remeshed_mat/'];
%%
base_pair_file = [data_dir 'test_pairs.txt'];

pair_list = {};
fid = fopen(base_pair_file); % open the file
while ~feof(fid) % feof(fid) is true when the file ends
      textLineEntry = fgetl(fid); % read one line
      [pair,comma] = split(textLineEntry,'	');
      pair_list{end+1} = pair;
end
fclose(fid); % close the file

NN_type = 'fast';
InitialGuess = 'normal_derivatives';

DS_num_eigs = 10;

radii_factor = 0.5;

weight_Orthonormality = 1;
weight_Proper = 1;
weight_Bijectivity = 1;

num_LB_eigs = 120;
ZO_start = 20;
ZO_step = 5;

load([data_dir,'Landmarks/lmk_list.mat']);
NUM_LMKS = length(lmk_list);

Steklov_settings = compute_steklov_settings(...
                NUM_LMKS, NN_type,InitialGuess,DS_num_eigs,radii_factor,...
                weight_Orthonormality,weight_Proper,weight_Bijectivity,...
                num_LB_eigs,ZO_start,ZO_step);
            

save_dir_root = sprintf('./output/');
if ~exist(save_dir_root,'dir')
    mkdir(save_dir_root);
end
%%
thresh = 0:0.001:0.50; % evaluation thresholds
r_ls = fliplr([0.25,0.5,0.75]);
for PAIR_IDX=1:length(pair_list)
    %
    pair = pair_list{PAIR_IDX};
    src_shape_name = pair{1};
    tar_shape_name = pair{2};
    %
    load([mat_dir,src_shape_name,'_R_0.75.mat']);
    Src=struct;
    Src.SHAPE = SHAPE;
    Src.landmarks = landmarks;
    Src.gamma = Gamma;
    Src_landmarks = Src.landmarks(lmk_list);
    %
    consistent_lmks = 1:length(Src.landmarks);% all landmarks are consistent in this dataset
    eval_landmarks = setdiff(consistent_lmks,lmk_list); % the landmarks on which to evaluate
    
    for r_idx=1:length(r_ls)
        r = r_ls(r_idx);
        %
        load([mat_dir,tar_shape_name,sprintf('_R_%.2f.mat',r)]);
        Tar = struct;
        Tar.SHAPE = SHAPE;
        Tar.landmarks = landmarks;
        Tar.gamma = Gamma;
        %
        Tar_landmarks = Tar.landmarks(lmk_list);
        %
        eval_Src = Src.landmarks(eval_landmarks);
        eval_Tar = Tar.landmarks(eval_landmarks);
        %
        save_dir = sprintf('%s/%03d_lmks/R_%.2f/',save_dir_root,NUM_LMKS,r);
        if ~exist(save_dir,'dir')
            mkdir(save_dir);
        end;
        %
        % -----------------------------------------------------
 
        fprintf('[%03d lmks R %i/%i] processing shape pair %s-%s (%i/%i)\n',NUM_LMKS,r_idx,length(r_ls),src_shape_name,tar_shape_name,PAIR_IDX,length(pair_list));

        save_dir_maps = [save_dir 'p2p/'];
        save_dir_err = [save_dir 'err/'];

        if ~exist(save_dir_maps,'dir')
            mkdir(save_dir_maps);
        end
        if ~exist(save_dir_err,'dir')
            mkdir(save_dir_err);
        end

        out_err_filename = [save_dir sprintf('%s_%s_err.mat',src_shape_name,tar_shape_name)];
        out_curve_filename = [save_dir sprintf('%s_%s_curve.mat',src_shape_name,tar_shape_name)];
        out_time_filename = [save_dir sprintf('%s_%s_time.mat',src_shape_name,tar_shape_name)];

        out_maps_filename = [save_dir_maps '%s_to_%s_p2p.mat'];
        out_err1_filename = [save_dir_err '%s_to_%s_err.mat'];
        out_curve1_filename = [save_dir_err '%s_to_%s_curve.mat'];

        if ~exist(out_err_filename,'file')||~exist(out_curve_filename,'file')||~exist(out_time_filename,'file')

            %% Adding the current radii to the Steklov settings struct
            Steklov_settings.landmarks_radii = computeLandmarkRadii(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings); %Radii of the disks at the landmarks.
            %% Compute Steklov elements
            tic;
            [Src_refined,Tar_refined,fullp2pTarSrc, ...
            fullp2pSrcTar,fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = ...
            compute_steklov(Src, Src_landmarks, ...
            Tar, Tar_landmarks, Steklov_settings);
            computation_time = toc/2.;
            %% Evaluation of the quality of the maps
            [err1,curve1,err2,curve2] = compute_bidirectional_err(Src,Tar,eval_Src,eval_Tar,fullp2pSrcTar_ZO,fullp2pTarSrc_ZO,thresh);
            curve_ZO = mean([curve1 ; curve2]);
            err_ZO = mean([err1 ; err2]);

            fprintf('\t %s to %s mean err = %.2e\n',src_shape_name,tar_shape_name,mean(err1));
            fprintf('\t %s to %s mean err = %.2e\n',tar_shape_name,src_shape_name,mean(err2));
            fprintf('\t %s-%s mean err = %.2e\n',src_shape_name,tar_shape_name,err_ZO);

            p2p = fullp2pSrcTar_ZO;
            err = err1;
            curve = curve1;
            save(sprintf(out_maps_filename,src_shape_name,tar_shape_name),'p2p');
            save(sprintf(out_err1_filename,src_shape_name,tar_shape_name),'err');
            save(sprintf(out_curve1_filename,src_shape_name,tar_shape_name),'curve');
            p2p = fullp2pTarSrc_ZO;
            err = err2;
            curve = curve2;
            save(sprintf(out_maps_filename,tar_shape_name,src_shape_name),'p2p');
            save(sprintf(out_err1_filename,tar_shape_name,src_shape_name),'err');
            save(sprintf(out_curve1_filename,tar_shape_name,src_shape_name),'curve');

            save(out_err_filename,'err_ZO');
            save(out_curve_filename,'curve_ZO');
            save(out_time_filename,'computation_time');
        end
    end
end