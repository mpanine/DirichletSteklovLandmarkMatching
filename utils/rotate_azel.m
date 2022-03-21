function [ VERT ] = rotate_azel( VERT,azel, alpha )
    %ROTATE_AZEL Summary of this function goes here
    %   Detailed explanation goes here
    % find unit vector for axis of rotation
    if numel(azel) == 2 % theta, phi
        theta = pi*azel(1)/180;
        phi = pi*azel(2)/180;
        u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
    elseif numel(azel) == 3 % direction vector
        u = azel(:)/norm(azel);
    end

    alph = alpha*pi/180;
    cosa = cos(alph);
    sina = sin(alph);
    vera = 1 - cosa;
    x = u(1);
    y = u(2);
    z = u(3);
    rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
           x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
           x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
       
    VERT = VERT*rot;
end

