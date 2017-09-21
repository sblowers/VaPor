disp('%%%%%% Establishing 1D-3D Domain Intersections %%%%')
disp('Establishing...')
[MdotArt,Vessel2Volume1,Volume2Vessel1,Mdot1] = vesselDomainExchange(Vessel1,L1,DomTot,GM_WM,BloodFlow,Option_BranchTerminationsOnly);
% Establish mass transfer between arterial domain segments and blood domain
% voxels.

[MdotVein,Vessel2Volume2,Volume2Vessel2,Mdot2] = vesselDomainExchange(Vessel2,L2,DomTot,GM_WM,-BloodFlow,Option_BranchTerminationsOnly);
% Establish mass transfer between venous domain segments and blood domain
% voxels.
if Option_CounterCurrentFlow
    Mdot = Mdot1-(Perfusion*VoxelSize^3);
    MdotB = (Perfusion*VoxelSize^3)+Mdot2;
else
    Mdot = Mdot1+Mdot2;
end
% Overall inter-domain mass transfer is the sum of both transfer from
% arteries and to veins.

disp('Complete.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
