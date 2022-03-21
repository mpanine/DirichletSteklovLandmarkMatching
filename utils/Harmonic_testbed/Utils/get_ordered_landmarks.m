function [ landmarks ] = get_ordered_landmarks( dataset, shape_name )
    %GET_LANDMARKS Returns the landmarks correctly ordered for each
    %dataset.
    lmks_fun = str2func(['order_' dataset '_lmks']);
    
    o = lmks_fun(shape_name);
    
    landmarks = o.landmarks;
end

