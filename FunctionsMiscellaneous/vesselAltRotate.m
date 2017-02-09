function VesselNew = vesselAltRotate(VesselOld,Axis,Angle)

VesselNew = VesselOld;

% find centroid
A = 0.5*(max(VesselOld(:,3)) - min(VesselOld(:,3))) + min(VesselOld(:,3));
B = 0.5*(max(VesselOld(:,4)) - min(VesselOld(:,4))) + min(VesselOld(:,4));
C = 0.5*(max(VesselOld(:,5)) - min(VesselOld(:,5))) + min(VesselOld(:,5));

% obtain axis
if strcmp(Axis,'x')
    U = 1;
    V = 0;
    W = 0;
elseif strcmp(Axis,'y')
    U = 0;
    V = 1;
    W = 0;
elseif strcmp(Axis,'z')
    U = 0;
    V = 0;
    W = 1;
else
    disp('error: axis argument must be string "x", "y", or "z"')
    return
end


L = U^2 + V^2 + W^2;
               
for n = 1:size(VesselOld,1)
        X = VesselOld(n,3);
        Y = VesselOld(n,4);
        Z = VesselOld(n,5);
        
        VesselNew(n,3) = ((A*(V^2+W^2)-U*(B*V+C*W-U*X-V*Y-W*Z))*(1-cosd(Angle))+L*X*cosd(Angle)+sqrt(L)*(-C*V+B*W-W*Y+V*Z)*sind(Angle))/L;
        VesselNew(n,4) = ((B*(U^2+W^2)-V*(A*U+C*W-U*X-V*Y-W*Z))*(1-cosd(Angle))+L*Y*cosd(Angle)+sqrt(L)*(C*U-A*W+W*X-U*Z)*sind(Angle))/L;
        VesselNew(n,5) = ((C*(U^2+V^2)-W*(A*U+B*V-U*X-V*Y-W*Z))*(1-cosd(Angle))+L*Z*cosd(Angle)+sqrt(L)*(-B*U+A*V-V*X+U*Y)*sind(Angle))/L;
end