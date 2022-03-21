function f = evalWave(S1, coeffs, decayrate, t)
%EVALWAVE evaluates the wave at time t
%S1: shape
%coeffs: coefficients of the wave at t=0 in the LB basis
%decayrate: rate of decay of the wave
%t: desired time

f = zeros(S1.nv,1);

for i = 1:length(coeffs)

    f = f + coeffs(i) * cos(S1.evals(i) * t) * exp(-decayrate*S1.evals(i)*t) * S1.evecs(:,i);

    
end




end

