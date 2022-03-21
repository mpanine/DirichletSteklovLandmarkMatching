function [ mean_val ] = meanError( entries, field )
    %MEANERROR Computes the mean error of the field for all entries.
    
    mean_val = 0;
    num_entries = length(entries);
    for i=1:num_entries
        mean_val = entries{i}.(field).mean/num_entries + mean_val;
    end
    
    
end

