function H = hash(Data)
%DIGEST Generates a hash based on the input struct
H = core.core_hash(Data);
H = core.b2h(H);   % To hex string
end

