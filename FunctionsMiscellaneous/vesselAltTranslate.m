function VesselNew = vesselAltTranslate(VesselOld,Translation)

VesselNew = VesselOld;

VesselNew(:,3) = VesselNew(:,3) + Translation(1);
VesselNew(:,4) = VesselNew(:,4) + Translation(2);
VesselNew(:,5) = VesselNew(:,5) + Translation(3);