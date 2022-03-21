classdef Settings
    %SETTINGS Base class for all Settings classes
    %Defines the hash update
    
    properties
        hash % The hexadecimal char array encoding the values of the properties within an instance of a Settings object.
    end
    
    methods
        function self = Settings()
            self.hash = '0';
        end
        
        function self = update(self)
            props = properties(self);
            
            hash = '';
            for i = 1:length(props)
              cur_prop = props{i};
              if ~strcmp(cur_prop, 'hash')
                value = self.(cur_prop);
                
                % Taking into account the property name for the hash
                % so that objects with same values don't xor themselves out
                hash_val = bitxor(core.h2b(core.hash(cur_prop),128),...
                                core.h2b(core.hash(value),128));
                
                % concatenating the hashes of each property
                hash = strcat(hash, core.b2h(hash_val));
              end
            end
            % Performing an MD5 digest on the concatenated hash
            % Reduces the probability of collisions
            hash = core.hash(hash);
            self.hash = hash;
        end
    end
end

