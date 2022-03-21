% Interactive landmark-based harmonic uniformization.
% Preserves the first three landmarks exactly. The rest are placed on the
% sphere automatically. -- THIS IS FALSE

%% Technical stuff

clear all;
close all;

addpath(genpath(pwd))

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings('a','b', false);


%% Load source
Src = struct;

Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/381.off', ...
                        shape_settings, cache_settings);



%% Load target
Tar = struct;

% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/388.off', ...
%                         shape_settings, cache_settings);

Tar.SHAPE = shape.compute('./Data/Faust_remeshed/vtx_5k/tr_reg_022.off', ...
                        shape_settings, cache_settings);

%%
global landmarks;
clear_landmarks();

%close all;
f1 = figure(1); 
set(f1,'Position',[100 100 200 400])

create_gui_Principled_Steklov(Src, Tar);

figure(2); 

hdt = datacursormode;
set(hdt,'UpdateFcn',{@read_landmarks,hdt,Src.SHAPE,Tar.SHAPE});
hold all;

subplot(1,2,1);
plot_function_faust(Tar.SHAPE.surface);
subplot(1,2,2);
plot_function_faust(Src.SHAPE.surface);
