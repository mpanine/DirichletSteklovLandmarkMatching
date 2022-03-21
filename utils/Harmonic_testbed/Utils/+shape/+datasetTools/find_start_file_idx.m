function start_file_idx = find_start_file_idx( files )
    continueSearch = true;
    start_file_idx = 3;
    while continueSearch
        if ~files(start_file_idx).isdir
            continueSearch = false;
        elseif start_file_idx > size(files, 1)
            continueSearch = false;
            throw(MException('GENERATELANDMARKFILES:directoryWithoutFiles',...
                'The data directory must contain at least one mesh file.'));
        else
            start_file_idx = start_file_idx + 1;
        end
    end
end

