function [ to_keep ] = random_exclude_3neighs( W, landmark_indices, NUM_LMKS )
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    random_indices = [];
    
    landmark_indices_orig = landmark_indices';
    
    for idx=1:NUM_LMKS
        % Random index
        selected_random = randperm(length(landmark_indices),1);
        landmark_index_random = landmark_indices(selected_random);

        % 2-ring neighs of indices
        neighs = find(W(landmark_index_random,:));
        landmark_indices = setdiff(landmark_indices,neighs);
        neighs2 = find(sum(W(neighs,:),1));
        landmark_indices = setdiff(landmark_indices,neighs2);
%         neighs2 = find(sum(W(neighs2,:),1));% 3-ring
%         % remove them from available indices
%         landmark_indices = setdiff(landmark_indices,neighs2);
        
        random_indices(end+1) = landmark_index_random;
    end
    
    if size(landmark_indices_orig,1)>size(landmark_indices_orig,2)
        landmark_indices_orig = landmark_indices_orig';
    end
    if size(random_indices,1)>size(random_indices,2)
        random_indices = random_indices';
    end
    
    to_keep = knnsearch(landmark_indices_orig',random_indices');
    
    % Replace indices too close with 3-neigh
    
end

