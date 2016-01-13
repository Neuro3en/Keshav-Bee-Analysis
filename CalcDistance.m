%% Distance formula
% calculates the euclidean distance between two coordinates in space
% Input:
% The cordinates both x and y for the two positions for which the distance
% is to be calculated

function euclideanDistance = CalcDistance(x1, y1, x2, y2)

euclideanDistance = sqrt((x2-x1)^2+(y2-y1)^2);
end

