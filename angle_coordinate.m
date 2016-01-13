%% function to find the angle between tow coordinates.
function ANGLE=angle_coordinate(x1,y1,x2,y2)
% Uses the equation of the staright line in the standard form, y=mx+c.
% where m is the slope and the tangent inverse of the slope gives the
% angle.
m=(y2-y1)/(x2-x1);
ANGLE=atan(m);

end
