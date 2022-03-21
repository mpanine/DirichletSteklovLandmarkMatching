function [R] = rot(axis,angle)
    R = eye(3);
    
    angle = deg2rad(angle);
    
    switch axis
        case 'x'
            R(2,:) = [0,cos(angle),-sin(angle)];
            R(3,:) = [0,sin(angle),cos(angle)];
        case 'y'
            R(1,:) = [cos(angle), 0, sin(angle)];
            R(3,:) = [-sin(angle),0,cos(angle)];
        case 'z'
            R(1,:) = [cos(angle),-sin(angle),0];
            R(2,:) = [sin(angle),cos(angle),0];
    end;
end