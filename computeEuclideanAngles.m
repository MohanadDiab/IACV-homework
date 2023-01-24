function angle = computeEuclideanAngles(l_a, l_b)
    %COMPUTEEUCLIDEANANGLES Returns the angle between two lines
    %   Given two lines, the function returns the angle between them
    angle = acos((l_a(1)*l_b(1)+l_a(2)*l_b(2))/(sqrt((l_a(1)^2+l_a(2)^2)*(l_b(1)^2+l_b(2)^2))));
    angle = 180*angle/pi;
end