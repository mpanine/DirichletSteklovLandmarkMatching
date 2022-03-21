% TESTS THE AREA REFINEMENT

clc; close all; clear all;
addpath(genpath('./'))


mesh_dir = 'data/';

k1 = 100;
k2 = k1;

meshOptions = {'IfComputeLB',false,'numEigs',k1,...
    'IfComputeGeoDist',false,... 
    'IfComputeNormals',true};

% s1_name = 'tr_reg_388';
s1_name = 'tr_reg_381';
% s2_name = 'tr_reg_388';
s1_index = 381; %Rotation correction indices. set to 0 for no rotation
s2_index = 388;

% s1_name = 'tr_reg_385'; % Giraffe and donkey
% s2_name = 'tr_reg_390';
% s1_index = 385; %Rotation correction indices. set to 0 for no rotation
% s2_index = 390;

% s1_name = 'tr_reg_386'; % Horse and coarse bull
% s2_name = 'tr_reg_397';
% s1_index = 386; %Rotation correction indices. set to 0 for no rotation
% s2_index = 397;


%% preprocess the mesh or load from the cache
S1 = MESH.preprocess([mesh_dir,s1_name], meshOptions{:});

preX = [S1.surface.X S1.surface.Y S1.surface.Z] * revert_SHREC07_rot( s1_index );
        
    S1.surface.X = preX(:,1);
    S1.surface.Y = preX(:,2);
    S1.surface.Z = preX(:,3);
    S1.surface.VERT  = preX;

%% Area refinement

area_factor = 3;
max_iterations = Inf;
S2 = faceAreaRefinement(S1, area_factor, max_iterations);


%% N-circle Neighborhood

num_circles = 2;
center = 5776;

[tri1, neigh1] = findTriangleNeigborhood(S1.surface, center, num_circles);

[tri2, neigh2] = findTriangleNeigborhood(S2.surface, center, num_circles);



%% Refine around point

S3 = S2;
S3.surface = refineVertexNeighborhood(S3.surface, center, num_circles);


%% Further refine around point

extra_ref = 5;

S4 = S3;

for r =1:extra_ref
    
   S4.surface = refineVertexNeighborhood(S4.surface, center, num_circles); 
    
end

%% Plotting
figure; 

subplot(1,4,1);
plot_function_faust(S1.surface);
hold on
plot3(S1.surface.VERT(neigh1,1),S1.surface.VERT(neigh1,2),S1.surface.VERT(neigh1,3),'.','color','r','markersize',20)
plot3(S1.surface.VERT(center,1),S1.surface.VERT(center,2),S1.surface.VERT(center,3),'.','color','b','markersize',25)
% plot3(S1.surface.VERT(1510,1),S1.surface.VERT(1510,2),S1.surface.VERT(1510,3),'.','color','r','markersize',20)
% plot3(S1.surface.VERT(5776,1),S1.surface.VERT(5776,2),S1.surface.VERT(5776,3),'.','color','r','markersize',20)
% plot3(S1.surface.VERT(1441,1),S1.surface.VERT(1441,2),S1.surface.VERT(1441,3),'.','color','r','markersize',20)
% plot3(S1.surface.VERT(1510,1),S1.surface.VERT(1510,2),S1.surface.VERT(1510,3),'.','color','r','markersize',20)
% plot3(S1.surface.VERT(5776,1),S1.surface.VERT(5776,2),S1.surface.VERT(5776,3),'.','color','r','markersize',20)
title('Original')

subplot(1,4,2);
plot_function_faust(S2.surface);
hold on
plot3(S2.surface.VERT(neigh2,1),S2.surface.VERT(neigh2,2),S2.surface.VERT(neigh2,3),'.','color','r','markersize',20)
plot3(S2.surface.VERT(center,1),S2.surface.VERT(center,2),S2.surface.VERT(center,3),'.','color','b','markersize',25)
% plot3(S2.surface.VERT(1441,1),S2.surface.VERT(1441,2),S2.surface.VERT(1441,3),'.','color','r','markersize',20)
% plot3(S2.surface.VERT(1510,1),S2.surface.VERT(1510,2),S2.surface.VERT(1510,3),'.','color','r','markersize',20)
% plot3(S2.surface.VERT(5776,1),S2.surface.VERT(5776,2),S2.surface.VERT(5776,3),'.','color','r','markersize',20)
title('Area Refined')

subplot(1,4,3);
plot_function_faust(S3.surface);
hold on
plot3(S3.surface.VERT(neigh2,1),S3.surface.VERT(neigh2,2),S3.surface.VERT(neigh2,3),'.','color','r','markersize',20)
plot3(S3.surface.VERT(center,1),S3.surface.VERT(center,2),S3.surface.VERT(center,3),'.','color','b','markersize',25)
% plot3(S2.surface.VERT(1441,1),S2.surface.VERT(1441,2),S2.surface.VERT(1441,3),'.','color','r','markersize',20)
% plot3(S2.surface.VERT(1510,1),S2.surface.VERT(1510,2),S2.surface.VERT(1510,3),'.','color','r','markersize',20)
% plot3(S2.surface.VERT(5776,1),S2.surface.VERT(5776,2),S2.surface.VERT(5776,3),'.','color','r','markersize',20)
title('Refined around Center Once')

subplot(1,4,4);
plot_function_faust(S4.surface);
hold on
plot3(S4.surface.VERT(neigh2,1),S4.surface.VERT(neigh2,2),S4.surface.VERT(neigh2,3),'.','color','r','markersize',20)
plot3(S4.surface.VERT(center,1),S4.surface.VERT(center,2),S4.surface.VERT(center,3),'.','color','b','markersize',25)
% plot3(S2.surface.VERT(1441,1),S2.surface.VERT(1441,2),S2.surface.VERT(1441,3),'.','color','r','markersize',20)
% plot3(S2.surface.VERT(1510,1),S2.surface.VERT(1510,2),S2.surface.VERT(1510,3),'.','color','r','markersize',20)
% plot3(S2.surface.VERT(5776,1),S2.surface.VERT(5776,2),S2.surface.VERT(5776,3),'.','color','r','markersize',20)
title(sprintf('Refined around Center %d times',extra_ref +1))

