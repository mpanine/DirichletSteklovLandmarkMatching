function E = centralFEenergy(S1, Basis2, landmarks1, landmarks2, NumControlLandmarks, radii1)
%CENTRALFEENERGY computes the control landmark-based energy


        S1 = insertCentralFElist(S1, landmarks1(1:end-NumControlLandmarks), radii1);
%         S2 = insertCentralFElist(S2, landmarks2(1:end-NumControlLandmarks), radii2);

        Basis1 = ComputeCentralHarmonicBasis(S1, landmarks1(1:end-NumControlLandmarks));
%         Basis2 = ComputeCentralHarmonicBasis(S2, landmarks2(1:end-NumControlLandmarks));


        E = norm( Basis1( landmarks1(end-NumControlLandmarks + 1:end), : ) - Basis2(landmarks2(end-NumControlLandmarks + 1:end), : ) )^2;




end