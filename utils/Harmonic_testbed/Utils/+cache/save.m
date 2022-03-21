function [ success ] = save( object, name, cache_settings )
%SAVE saves the object according to the cache settings.
%   Inputs:
%       - object: the object to cache
%       - name: the name given to the cached object
%       - cache_settings: the settings of the cache
%   Output: success (logical)
    
    if cache_settings.enabled
        if isdir(cache_settings.directory)
            if cache_settings.verbose
                fprintf(sprintf('Saving to cache object with settings %s...', class(cache_settings))); tic;
            end
            saveDir = sprintf('%s/%s/%s/', cache_settings.directory, ...
                                           class(object.settings),...
                                           object.settings.hash);
            if ~exist(saveDir, 'dir')
                mkdir(saveDir)
                settingsInCache = object.settings;
                save(sprintf('%s/%s.mat',saveDir, 'settings'),'settingsInCache');
                clear settingsInCache;
            end
            save(sprintf('%s/%s.mat',saveDir, name),'object');
            if cache_settings.verbose
                t = toc; fprintf('done: %.4f\n', t);
            end
            success = true;
        else
            success = false;
        end
    else
        success = false;
        if cache_settings.verbose
            fprintf('Cache is disabled. Set cache_settings.enabled to true.\n');
        end
    end
end

