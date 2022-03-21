function [ S ] = revert_FAUST_rot( S )
%CORRECT_FAUST_ROTATION rectifies the rotation of FAUST models so that each
%character is upright.
    S.surface.VERT = S.surface.VERT*rot('x',-90);
    S.surface.X = S.surface.VERT(:,1);
    S.surface.Y = S.surface.VERT(:,2);
    S.surface.Z = S.surface.VERT(:,3);
end

