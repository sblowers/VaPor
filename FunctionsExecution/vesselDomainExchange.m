function [MdotVessel,Vessel2Volume,Volume2Vessel,MdotVoxels] = vesselDomainExchange(Vessel,L,DomainFull,DomainSource,BloodFlow)


%%%%%% Create Inter-Domain Cells %%%%%%%%%%%%%%%%
Vessel2Volume = cell(size(Vessel,1),1);
Volume2Vessel = cell(size(DomainFull));
% Create cells that store information about inter-domain transfer. For each
% voxel, Volume2Vessel contains information about intersecting vessel
% segments for each voxel and how much transfer occurs between them.
% Vessel2Volume contains information about intersecting voxels for each
% vessel segment and how much transfer occurs between them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initial Guess For Inter-Domain Transfer %%
MdotVessel = BloodFlow * L/sum(L);
% Inter-Domain transfer based on proportional length of segment. This is
% subject to alteration further on if certain segments do not have voxel
% intersections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Find Intersections for Segments %%%%%%%%%%
for iv = 2:size(Vessel,1)
    Intersections = cubeIntersect(Vessel(iv,3:5),Vessel(Vessel(iv,7),3:5));
    % Find intersections between vessels and cubes.
    
    if ~isempty(Intersections)
        IntersectionsBool = true(size(Intersections,1),1);
        for iw = 1:size(Intersections,1)
            if ~DomainFull(Intersections(iw,1),Intersections(iw,2),Intersections(iw,3))
                IntersectionsBool(iw) = false;
            end
        end
        Intersections = Intersections(IntersectionsBool,:);
    end
    % Deleting intersections not within domain.
    
    
    if isempty(Intersections) && (round(Vessel(iv,3))>=1 && round(Vessel(iv,3))<=size(DomainFull,1))...
            && (round(Vessel(iv,4))>=1 && round(Vessel(iv,4))<=size(DomainFull,2))...
            && (round(Vessel(iv,5))>=1 && round(Vessel(iv,5))<=size(DomainFull,3))
        if DomainFull(round(Vessel(iv,3)),round(Vessel(iv,4)),round(Vessel(iv,5)))
            Intersections = [round(Vessel(iv,3)),round(Vessel(iv,4)),round(Vessel(iv,5))];
        end
    end
    % Check if vessel is within a cube (but doesn't intersect the boundries
    % of that cube).
    
    if ~isempty(Intersections)
        Intersections(:,end+1) = 0;
    end
    % Pre-allocate source term to any intersections found.
    
    Vessel2Volume{iv} = Intersections;
    % Update cell.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Voxel Inter-Domain Transfer %%%%%%%%%%%%%%
MdotVoxels = zeros(size(DomainFull));
% Preallocate array for inter-domain mass transfer in voxels.

for iv = 2:size(Vessel,1) % For every line segment.
    if ~isempty(Vessel2Volume{iv}) % If there is an intersection.
        for  iw=1:size(Vessel2Volume{iv},1)
            i = Vessel2Volume{iv}(iw,1);
            j = Vessel2Volume{iv}(iw,2);
            k = Vessel2Volume{iv}(iw,3); % Get location of intersection voxel.
            if DomainSource(i,j,k) % If intersection voxel is in domain that allows inter-domain transfer.
                Vessel2Volume{iv}(iw,4) = MdotVessel(iv)/numel(Vessel2Volume{iv}(iw,:)); % Update Vessel2Volume to include this specific transfer.
                MdotVoxels(i,j,k) = MdotVoxels(i,j,k) + MdotVessel(iv)/numel(Vessel2Volume{iv}(iw,:)); % Update total transfer to voxel.
            end
        end
    end
end
% This allocates the transfers from each segement to the corresponding
% voxels. If there are any segments that do not transfer then the result
% will be unbalanced. This is then corrected in the subsequent section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Adjust Inter-Domain Transfer %%%%%%%%%%%%%
if sum(sum(sum(MdotVoxels))) == 0
    error('Error: No intersections between vessel tree and domain found.')
    % Check if there is any transfer to any voxel. Otherwise domains are
    % not connected and program cannot proceed.
else
    FracSource = BloodFlow/sum(sum(sum(MdotVoxels))); % Amount MdotVessels need to be adjusted to balance the inter-domain transfer.
    MdotVoxels = MdotVoxels * FracSource; % Adjusting the voxel inter-domain transfer to match overall BloodFlow.
    for iv = 1:size(Vessel,1) % For every line segment.
        if ~isempty(Vessel2Volume{iv}) % If there is an intersection.
            Vessel2Volume{iv}(:,4) = Vessel2Volume{iv}(:,4)*FracSource; % Adjusting the segment inter-domain transfer so that overall mass transfer is conserved.
            MdotVessel(iv) = sum(Vessel2Volume{iv}(:,4)); % Overall segment transfer is sum of flow to each voxel.
            for iw = 1:size(Vessel2Volume{iv},1) % For every voxel that this segment intersects.
                ijk = Vessel2Volume{iv}(iw,1:3);
                Volume2Vessel{ijk(1),ijk(2),ijk(3)}(end+1,:) = [iv,Vessel2Volume{iv}(iw,4)]; % Update stored value for inter-domain transfer in Volume2Vessel.
            end
        else
            MdotVessel(iv) = 0; % If there are no intersections segment inter-domain flow is zero.
        end
    end
end
% As some vessels do not transfer blood to their corresponding
% intersections, the overall transfer becomes unbalanced. This part fixes
% the imbalance by setting transfer in appropriate segments to zero and
% proportionally increasing inter-domain flow in the segments that do
% transfer flow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end