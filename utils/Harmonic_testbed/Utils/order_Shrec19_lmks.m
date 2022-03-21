function [ landmarks_obj ] = order_Shrec19_lmks( shape_name )
    %ORDER_SHREC19_LMKS Summary of this function goes here
    %   Detailed explanation goes here
    
%     o = load(['./Data/Shrec19/Landmarks_FARM_first/' shape_name]);
%     first_lmks = o.landmarks;
%     
%     landmarks_obj = load(['./Data/Shrec19/Landmarks/' shape_name]);
%     first_lmks_idx = knnsearch(landmarks_obj.landmarks,first_lmks);
%     landmarks_obj.landmarks(first_lmks_idx) = [];
%     landmarks_obj.landmarks = [first_lmks; landmarks_obj.landmarks];

    landmarks_obj = load(['./Data/Shrec19/Landmarks_Re-Ordered/' shape_name]);
    
end

