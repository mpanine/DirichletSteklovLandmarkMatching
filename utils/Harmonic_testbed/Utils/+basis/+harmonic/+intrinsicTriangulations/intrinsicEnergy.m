function E = intrinsicEnergy(shape_number, S1, Basis2, landmarks1, landmarks2, NumControlLandmarks, landmarks_radii)
%CENTRALFEENERGY computes the control landmark-based energy



        Basis1 = computeIntrinsicBasis( shape_number, S1.nv, landmarks1(1:end-NumControlLandmarks), landmarks_radii );

        E = norm( Basis1( landmarks1(end-NumControlLandmarks + 1:end), : ) - Basis2(landmarks2(end-NumControlLandmarks + 1:end), : ) )^2;




end