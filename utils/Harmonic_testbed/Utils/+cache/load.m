function [ cached_data ] = load( cache_settings, name, settings )
%LOAD Retrieves the cache based on the cache directory and the settings
%struct associated to the cached data.
%   Inputs:
%   - cache_settings: settings of the cache
%   - name: name of the shape
%   - settings: settings struct of the cached data
%   Output: cached_data
    if cache_settings.enabled
        % Generate a hash string corresponding to the settings
        hashDir = sprintf('%s/%s/', class(settings), settings.hash);
        save_filename = [cache_settings.directory,hashDir,name,'.mat'];
        if exist(save_filename,'file')
            if cache_settings.verbose
                fprintf(sprintf('Loading from cache object with settings %s...', class(settings))); tic;
            end
            cached_data = load(save_filename);
            if cache_settings.verbose
                t = toc; fprintf('done: %.4fs\n',t);
            end
        else
            cached_data = struct;
            if cache_settings.verbose
                fprintf('Could not load cache for shape %s.\n',name);
            end
        end
    else
        cached_data = struct;
        if cache_settings.verbose
            fprintf('Cache is disabled. Set cache_settings.enabled to true.\n');
        end
    end
end

