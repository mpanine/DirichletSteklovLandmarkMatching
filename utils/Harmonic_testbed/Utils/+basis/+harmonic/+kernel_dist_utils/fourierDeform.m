function defB = fourierDeform(Basis, coeffs)
%FOURIERDEFORM deforms the basis w.r.t Fourier coefficients
%Basis: basis of harmonic functions
%coeffs: Fourier coefficients of the deformation size = (num_Basis_funcs, highest_fourier_order)


defB = Basis;

for i = 1:size(coeffs,2)

    defB = defB + repmat(coeffs(:,i)', [size(Basis, 1) 1]).* sin(pi * i * Basis);


end




end

