function [X,Y,Z] = createCylinder(Pnt1,Pnt2,R1,R2)

% Create a cylinder of 10 points in length
if abs(R1-R2) <= 1e-5
    [X,Y,Z] = cylinder(ones(10,1)*R1,10); 
    % If cylinder diameter difference between points is less than 0.01mm 
    % then assume they are equal diameters.
else
    [X,Y,Z] = cylinder(R1:(R2-R1)/9:R2,10);
    % If cylinder diameter difference between points is less than 0.01mm 
    % then assume they are equal diameters.
end


% Elongate
Z = Z*norm(Pnt1-Pnt2);

% Translate
X = X + Pnt1(1);
Y = Y + Pnt1(2);
Z = Z + Pnt1(3);

% Rotate

Angle = acosd(dot(Pnt2-Pnt1,[0,0,1])/(norm(Pnt2-Pnt1)*(norm([0,0,1]))));
CrossPro = cross([0,0,1],Pnt2-Pnt1);

A = Pnt1(1);
B = Pnt1(2);
C = Pnt1(3);
U = CrossPro(1);
V = CrossPro(2);
W = CrossPro(3);
L = U^2 + V^2 + W^2;

X2 = X;
Y2 = Y;
Z2 = Z;

if abs(Angle) > 1e-5
    if abs(U)<1e-5 && abs(V)<1e-5 && abs(W)<1e-5
        Z = Pnt1(3)-(Z-Pnt1(3)); % for rotating 180 degrees
    else
        X = ((A*(V^2+W^2)-U*(B*V+C*W-U*X2-V*Y2-W*Z2)).*(1-cosd(Angle))+L*X2*cosd(Angle)+sqrt(L).*(-C*V+B*W-W*Y2+V*Z2)*sind(Angle))/L;
        Y = ((B*(U^2+W^2)-V*(A*U+C*W-U*X2-V*Y2-W*Z2)).*(1-cosd(Angle))+L*Y2*cosd(Angle)+sqrt(L).*(C*U-A*W+W*X2-U*Z2)*sind(Angle))/L;
        Z = ((C*(U^2+V^2)-W*(A*U+B*V-U*X2-V*Y2-W*Z2)).*(1-cosd(Angle))+L*Z2*cosd(Angle)+sqrt(L).*(-B*U+A*V-V*X2+U*Y2)*sind(Angle))/L;
    end
else
    
end

end

