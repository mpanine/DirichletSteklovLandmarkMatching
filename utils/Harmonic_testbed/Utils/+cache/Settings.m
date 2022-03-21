classdef Settings < core.Settings
    %SETTINGS Settings for a collection of cached files
    
    properties
        name % The name of the cache.
        root % The path to the directory where to create the cache.
        directory % The path to the cache directory. It is equal to [***root***, ***name***, '/']
        enabled % If true, the caching mechanism is active.
        verbose % If true, all execution information texts are displayed.
    end
    
    methods
        function settings = Settings(enabled, name, root, verbose)
            if nargin<1
                enabled = false;
            end
            if nargin<2
                name = 'none';
            end
            if nargin<3
                root = sprintf('%s/Cache/',pwd);
            end
            if nargin<4
                verbose = false;
            end
            assert(ischar(name),'name has to be a string.');
            
            if enabled
                if ~exist(root, 'dir')
                   mkdir(root)
                end
                if ~exist(root, 'dir')
                   mkdir(root)
                   warning('Creating root directory for the cache at %s',root);
                end
            end
            
            settings.name = name;
            settings.root = root;
            settings.enabled = enabled;
            settings.verbose = verbose;
            
            settings.directory = sprintf('%s%s/', settings.root, settings.name);
            
            if enabled
                if ~exist(settings.directory, 'dir')
                   mkdir(settings.directory)
                end
            end
        end
        
        function self = set.name(self, value)
            assert(ischar(value), 'name has to be a string.');
            self.name = value;
            self = self.update();
        end
        
        function self = set.root(self, value)
            assert(ischar(value), 'root has to be a string.');
            self.root = value;
            self = self.update();
        end
        
        function self = set.enabled(self, value)
            assert(islogical(value), 'enabled has to be a logical.');
            self.enabled = value;
            self = self.update();
        end
    end
    
end

