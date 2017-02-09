function VesselNew = vesselAltScale(VesselOld,Scale,Axis)

if nargin < 3
    ScaleX = Scale;
    ScaleY = Scale;
    ScaleZ = Scale;
else
    if strcmp(Axis,'x')
        ScaleX = Scale;
        ScaleY = 1;
        ScaleZ = 1;
    elseif strcmp(Axis,'y')
        ScaleX = 1;
        ScaleY = Scale;
        ScaleZ = 1;
    elseif strcmp(Axis,'z')
        ScaleX = 1;
        ScaleY = 1;
        ScaleZ = Scale;
    else
        disp('error: axis argument must be string "x", "y", or "z"')
        return
    end
end

VesselNew = VesselOld;

A = 0.5*(max(VesselOld(:,3)) - min(VesselOld(:,3))) + min(VesselOld(:,3));
B = 0.5*(max(VesselOld(:,4)) - min(VesselOld(:,4))) + min(VesselOld(:,4));
C = 0.5*(max(VesselOld(:,5)) - min(VesselOld(:,5))) + min(VesselOld(:,5));

VesselNew(:,3) = VesselNew(:,3) - A;
VesselNew(:,4) = VesselNew(:,4) - B;
VesselNew(:,5) = VesselNew(:,5) - C;

VesselNew(:,3) = VesselNew(:,3) * ScaleX;
VesselNew(:,4) = VesselNew(:,4) * ScaleY;
VesselNew(:,5) = VesselNew(:,5) * ScaleZ;

VesselNew(:,3) = VesselNew(:,3) + A;
VesselNew(:,4) = VesselNew(:,4) + B;
VesselNew(:,5) = VesselNew(:,5) + C;
