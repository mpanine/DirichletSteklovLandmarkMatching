function Out = powersOfVector(v, max_pow)
%POWERSOFVECTOR builds a matrix of powers of a vector (1 to max_pow)
%v: vetor
%max_pow: highest desired power.

Out = zeros(length(v), max_pow);



for p = 1:max_pow
    
    Out(:, p) = v.^p;

end    


end

