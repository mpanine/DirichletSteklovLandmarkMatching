function [ landmarks_obj ] = order_Shrec07_Fourleg_lmks( shape_name )
    %ORDER_SHREC07_FOURLEG_LMKS 
    
%     first_lmks_idx = [14,10,8,12,6,13,15];
    
    landmarks_obj = load(['./Data/Shrec07_Fourleg/Landmarks_Re-Ordered/' shape_name]);
%     first_lmks = landmarks_obj.landmarks(first_lmks_idx);
%     landmarks_obj.landmarks(first_lmks_idx) = [];
%     landmarks_obj.landmarks = [first_lmks; landmarks_obj.landmarks];
    
end

