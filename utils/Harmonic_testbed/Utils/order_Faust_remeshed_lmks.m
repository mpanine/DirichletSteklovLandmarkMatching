function [ landmarks_obj ] = order_Faust_remeshed_lmks( shape_name )
    %ORDER_FAST_REMESHED_LMKS Summary of this function goes here
    %   Detailed explanation goes here
    
%     first_lmks_vertex_idx = [4254,4630,2135,3642,1149,114,2761,4109]';
%     first_lmks_idx = [1814,2659,783,3299,3141,2881,370,471]';
    first_lmks_idx = [2351,639,3284,3082,2709,2984,1950,850]';
    
    landmarks_obj = load(['./Data/Faust_remeshed/Landmarks/' shape_name]);
%     knnsearch(landmarks_obj.landmarks,first_lmks_idx)
    first_lmks = landmarks_obj.landmarks(first_lmks_idx);
    landmarks_obj.landmarks(first_lmks_idx) = [];
    landmarks_obj.landmarks = [first_lmks; landmarks_obj.landmarks];
    
    
end

