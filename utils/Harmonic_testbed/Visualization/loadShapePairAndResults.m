% Loads a pair of shapes from directory shapedir, as guided by the number
% of the shape pair pairnum and the results stored in directory resultdir.

clear all;

%% Settings


shapedir = './Data/TOSCA_Isometric/vtx_5k/';
resultdir = './Results/flexible_largescale/TOSCA_Isometric_MHB_ISO_ZO_6_LANDMARKS_PRE_REF/';
pairnum = 25;

load([resultdir 'settings.mat']);




%% Parsing the results to get the shapes to load

S = dir([resultdir '/' num2str(pairnum) '_*_C.mat']); %We use the Fmap file to parse for the shape names.
prenames = S.name(numel(num2str(pairnum))+2:end-6); %Remove the non-informative characters

names = split(prenames,'-'); % splits the names


%% Load shapes

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', true,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', true);
cache_settings = cache.Settings(dataset,'./Cache/');


Src.SHAPE = shape.compute([shapedir,'/',names{1},'.off'], ...
                        shape_settings, cache_settings);             
Tar.SHAPE = shape.compute([shapedir,'/',names{2},'.off'], ...
                        shape_settings, cache_settings);


%% Load Settings and Landmarks


landmarks_obj = load(['./Data/' dataset '/Landmarks/' names{1}]);%1:min(Src.SHAPE.nv,Tar.SHAPE.nv);%shape.getDatasetLandmarks( dataset, mesh_dir, shape_name_Src );     
landmarks_Src = landmarks_obj.landmarks;
landmarks_obj = load(['./Data/' dataset '/Landmarks/' names{2}]);%1:min(Src.SHAPE.nv,Tar.SHAPE.nv);%shape.getDatasetLandmarks( dataset, mesh_dir, shape_name_Tar );      
landmarks_Tar = landmarks_obj.landmarks;
                    
                    
                    
%% Load Maps


SS = dir([resultdir '/' num2str(pairnum) '_*.mat']);


for i = 1:length(SS)
    
    load([resultdir SS(i).name])
    
end




fprintf('\nPair Loaded.\n')




