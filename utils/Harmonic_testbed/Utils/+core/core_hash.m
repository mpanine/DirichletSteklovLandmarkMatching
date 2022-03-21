function [ H ] = core_hash( Data )
%CORE_HASH generates the core hash based on input data
H = core.h2b('0',128);
% Consider the type of empty arrays:
if isa(Data, 'struct')
    n = numel(Data);
    if n == 1  % Scalar struct:
        F = sort(fieldnames(Data));
        for iField = 1:length(F)
            H = bitxor(H, core.core_hash(F{iField})); % Taking into account the field's name
            H = bitxor(H, core.core_hash(Data.(F{iField})));
        end
    else  % Struct array:
      for iS = 1:n
         H = bitxor(H, core.core_hash(Data(iS)));
      end
   end
elseif isempty(Data)
   % No further actions needed
elseif isnumeric(Data)
    H = bitxor(H, core.h2b(core.GetMD5(full(Data)),128));
elseif ischar(Data)
    H = bitxor(H, core.h2b(core.GetMD5(Data),128));
elseif iscell(Data)
    for iS = 1:numel(Data)
       H = bitxor(H, core.core_hash(Data{iS}));
    end
elseif islogical(Data)
    H = bitxor(H, core.h2b(core.GetMD5(Data),128));
elseif isa(Data, 'function_handle')
    % NEW version (more lightweight):
    fun_struct = functions(Data);
    H = core.h2b(core.GetMD5(strcat(core.b2h(H),fun_struct.function)),128);
%     % OLD version:
%     H = bitxor(H, core.core_hash(functions(Data)));
elseif isprop(Data, 'hash')
    H = core.h2b(core.GetMD5(strcat(core.b2h(H),Data.hash)),128);
%     H = bitxor(H, core.h2b(Data.hash,128));
else
   warning(['Type of variable not considered: ', class(Data)]);
end

end
