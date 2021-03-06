function [ S ] = revert_SHREC07_rot( S )
%REVERT_SHREC07_ROT reverts the rotation applied to the models in the
%SHREC07 dataset

    index = str2num(S.name);
    
    switch index
        case 381
            Rotation_matrix = rot('y',90)*rot('x',-90);
        case 382
            Rotation_matrix = rot('z',180);
        case 383
            Rotation_matrix = eye(3);
        case 384
            Rotation_matrix = rot('x',-180);
        case 385
            Rotation_matrix = eye(3);
        case 386
            Rotation_matrix = rot('z',-90);
        case 387
            Rotation_matrix = rot('x',-75)*rot('z',-90);
        case 388
            Rotation_matrix = rot('x',-90)*rot('z',90);
        case 389
            Rotation_matrix = eye(3);
        case 390
            Rotation_matrix = rot('x',-90);
        case 391
            Rotation_matrix = rot('x',-90);
        case 392
            Rotation_matrix = rot('x',-90)*rot('z',90);
        case 393
            Rotation_matrix = eye(3);
        case 394
            Rotation_matrix = rot('z',90);
        case 395
            Rotation_matrix = rot('x',-90)*rot('z',-90);
        case 396
            Rotation_matrix = eye(3);
        case 397
            Rotation_matrix = rot('x',-90);
        case 398
            Rotation_matrix = rot('y',90)*rot('x',-90);
        case 399
            Rotation_matrix = rot('x',-90);
        case 400
            Rotation_matrix = rot('x',-90);
    end;
    
    S.surface.VERT = S.surface.VERT*Rotation_matrix;
    S.surface.X = S.surface.VERT(:,1);
    S.surface.Y = S.surface.VERT(:,2);
    S.surface.Z = S.surface.VERT(:,3);
end

