function [D2,T2] = distPointToLineSeg(Vec1,Vec2,Pnt)
% Finds shortest distance between a point (pnt) and a line defined
% by two points (vec1 vec2).
l_squared = sum((Vec2 - Vec1).^2);
if l_squared == 0
    D2 = norm(Pnt - Vec1);
else
    T2 = max(0, min(1, dot(Pnt - Vec1, Vec2 - Vec1) / l_squared));
    projection = Vec1 + T2 * (Vec2 - Vec1);
    D2 = norm(Pnt-projection);
end
end
