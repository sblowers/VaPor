
%%%%%%% Loading Vessels from Files %%%%%%%%%%%%%%
if Option_LoadPrevious
    disp('%%%%%% Loading Previous Results %%%%%%%%%%%%%%%%%%%')
    disp(['Loading: ' LoadPreviousFilename])
    load(LoadPreviousFilename,'Vessel1','Vessel2','InletPoints','OutletPoints')
    disp('Loaded')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    % If Options_LoadPrevious is used then the alterations to the vessel
    % tree do not need to be performed.
    
    if Option_Stroke
        StrokeLocation = findPoints(StrokeLocationPos,Vessel1);
    end
    % If Vessels are loaded, then this finds the correct location for
    % simulating the stroke.
    
else
    if strcmp(Option_VesselAdjustment, 'default')
        % Commands for extracting files
        Vessel1 = load([ArteriesLocation ArteriesFile]); % Load Arteries.
        
        Vessel2 = load([VeinsLocation VeinsFile],'Vessel'); % Load Veins.
        Vessel2 = Vessel2.Vessel; % Vessel2 is loaded as a structure.
        
        Vessel1 = vesselAltRotate(Vessel1,'x',90);
        Vessel1 = vesselAltScale(Vessel1,1.2,'z');
        Vessel1 = vesselAltScale(Vessel1,0.65);
        Vessel1 = vesselAltTranslate(Vessel1,[-25,140,-40]);
        % Alterations for the arterial tree to fit into domain.
        
        Vessel2 = vesselAltScale(Vessel2,256/100,'z'); % only 100 pixels in the z axis and 256 pixels in x and y.
        Vessel2 = vesselAltRotate(Vessel2,'x',60);
        Vessel2 = vesselAltScale(Vessel2,0.60);
        Vessel2 = vesselAltTranslate(Vessel2,[-70,-60,-5]);
        % Alterations for the venous tree to fit into domain.
    else
        Vessel1 = load([ArteriesLocation ArteriesFile]);
        Vessel2 = load([VeinsLocation VeinsFile]);
        % If another arterial or venous tree is used then any alterations
        % should be included here to fit the tree to the domain.
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Saving Original Diameters %%%%%%%%%%%%%%%%
Vessel1OriginalDiameters = Vessel1(:,6);
Vessel2OriginalDiameters = Vessel2(:,6);
% As the diameters are altered within the 1D flowrate derivations, these
% are saved for comparison if desired. If the vessel tree is altered
% through generation then the values might not match up to the original
% nodes in the new tree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Adjusting Vessels if Required %%%%%%%%%%%%
if Option_LoadPrevious
    if strcmp(Option_LoadPreviousCoarsening, 'two') || strcmp(Option_LoadPreviousCoarsening, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1)./2 + 1; % Scale arterial tree down by a factor of 2.
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1)./2 + 1; % Scale venous tree down by a factor of 2.
    end
    
    if strcmp(Option_LoadPreviousCoarsening, 'three') || strcmp(Option_LoadPreviousCoarsening, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1)./3 + 1; % Scale arterial tree down by a factor of 3.
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1)./3 + 1; % Scale venous tree down by a factor of 3.
    end
    % If Option_LoadPrevious is used then the corresponding coarsening is
    % defined by Option_LoadPreviousCoarsening and not by the domain
    % Option_DomainCoarsening.
    
    if strcmp(Option_LoadPreviousRefining, 'two') || strcmp(Option_LoadPreviousRefining, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1).*2 + 1; % Scale arterial tree down by a factor of 2.
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1).*2 + 1; % Scale venous tree down by a factor of 2.
    end
    
    if strcmp(Option_LoadPreviousRefining, 'three') || strcmp(Option_LoadPreviousRefining, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1).*3 + 1; % Scale arterial tree down by a factor of 3.
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1).*3 + 1; % Scale venous tree down by a factor of 3.
    end
    % If Option_LoadPrevious is used then the corresponding refining is
    % defined by Option_LoadPreviousRefining.
    
else
    
    if strcmp(Option_DomainCoarsening, 'two') || strcmp(Option_DomainCoarsening, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1)./2 + 1; % Scale arterial tree down by a factor of 2
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1)./2 + 1; % Scale venous tree down by a factor of 2
    end
    
    if strcmp(Option_DomainCoarsening, 'three') || strcmp(Option_DomainCoarsening, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1)./3 + 1; % Scale arterial tree down by a factor of 3
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1)./3 + 1; % Scale venous tree down by a factor of 3
    end
    % If the domain has been coarsened using Option_DomainCoarsening, then
    % the vessels are scaled to match the new voxel numbers.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Overwriting for a Cube Brain %%%%%%%%%%%%%
if Option_CubeBrain
    % If the Option_CubeBrain is used then the vessels are overwritten to
    % match the corresponding option.
    
    if Option_CubeBrainLinear
        Vessel1 = [1,0,0,2,2,0.01,0; 2,0,1,2,2,0.01,1; 3,0,2,2,2,0.01,2];
        Vessel2 = [1,0,0,2,CubeSize+1,0.01,0; 2,0,1,2,CubeSize+1,0.01,1; 3,0,2,2,CubeSize+1,0.01,2];
        % Creates vessels that intersect either end of the linear cube brain.
        
        InletPoints = 1;
        OutletPoints = 1;
        InletFlows = 1;
        OutletFlows = 1;
        % Set inlet points and flow weights for vasculature.
        
    elseif Option_CubeBrainArm || Option_CubeBrainArm2
        if Option_CubeBrainArm2
            Arm_Midpoint = 1.5;
        else
            Arm_Midpoint = round((CubeSize+2)/2)+0.5;
        end
        % Find and set central voxels of the arm.
        
        Arm_VesselEnds = round(CubeArmRatio*CubeSize)+0.5; % Find last voxel of arm.
        Vessel1 = zeros(Arm_NumVessels+1,7);
        Vessel1(:,1) = (1:Arm_NumVessels+1)';
        Vessel1(:,3) = Arm_Midpoint*ones(Arm_NumVessels+1,1);
        Vessel1(:,4) = Arm_Midpoint*ones(Arm_NumVessels+1,1);
        Vessel1(1:end-1,5) = linspace(0,Arm_VesselEnds,Arm_NumVessels)'; % Creates the number of vessel segments as specified by Arm_NumVessels.
        Vessel1(end,5) = Arm_VesselEnds + 1/10; % Ensures last segment is in last voxel of arm.
        Vessel1(:,6) = 0.01*ones(Arm_NumVessels+1,1);
        Vessel1(:,7) = (0:Arm_NumVessels)';
        % Creates arterial vessel that passes through the centre of the cubic
        % arm model created.
        
        if Option_CubeBrainArmCounterCurrent
            Vessel2 = Vessel1;
            % For counter-current arm. Vessels occupy the same central voxels.
        else
            Vessel2 = [1,0,Arm_Midpoint,Arm_Midpoint,Arm_VesselEnds+2,0.01,0;...
                2,0,Arm_Midpoint,Arm_Midpoint,Arm_VesselEnds+1,0.01,1;...
                3,0,Arm_Midpoint,Arm_Midpoint,Arm_VesselEnds+1/10,0.01,2;...
                4,0,Arm_Midpoint,Arm_Midpoint,Arm_VesselEnds,0.01,3];
            % The vein is an extension of the arteries that extends in the same
            % direction from where the arteries terminates along the central
            % voxels. Non-counter-current arm.
        end
        
        InletPoints = 1;
        OutletPoints = 1;
        InletFlows = 1;
        OutletFlows = 1;
        % Set inlet points and flow weights for vasculature.
        
    else
        Vessel1 = zeros(3,7);
        Vessel2 = zeros(3,7);
        Vessel1(:,1) = 1:3; Vessel1(:,2) = 1; Vessel1(:,3) = 2.5; Vessel1(:,4) = 2.5;
        Vessel1(:,5) = linspace(0,2.5,3); Vessel1(:,6) = 0.01; Vessel1(:,7) = 0:2;
        Vessel2(:,1) = 1:3; Vessel2(:,2) = 1; Vessel2(:,3) = CubeSize+0.5; Vessel2(:,4) = CubeSize+0.5;
        Vessel2(:,5) = linspace(CubeSize+1.5,CubeSize,3); Vessel2(:,6) = 0.01; Vessel2(:,7) = 0:2;
        % Creates an artery and a vein at opposite corners of the cube, three
        % segments in length.
        
        InletPoints = 1;
        OutletPoints = 1;
        InletFlows = 1;
        OutletFlows = 1;
        % Set inlet points and flow weights for vasculature.
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Generating Vessels if Required %%%%%%%%%%%
if Option_GenerateVessels
    InletOriginal = Vessel1(InletPoints,:);
    OutletOriginal = Vessel2(OutletPoints,:);
    % Storing the inlets and outlets with their locations so they can be
    % found later.
    
    if Option_Stroke
        StrokeOriginal = Vessel1(StrokeLocation,:);
    end
    % Storing the location of stroke (if used) so it can be found
    % afterwards.
    
    disp('%%%%%%% Generating Arterial Vessels %%%%%%%%%%%%%%%') % Display.
    tic % Start timing artery generation.
    if Option_RRTstar;
        Vessel1 = vesselGenerationRRT_star(Vessel1,Perfusion,GM_WM,ArteryGenerationIterations,ArteryGenerationWeightFactor,RRTstarEpsilon/VoxelSize); % Generating Arteries
    else
        Vessel1 = vesselGenerationRRT(Vessel1,Perfusion,GM_WM,ArteryGenerationIterations,ArteryGenerationWeightFactor); % Generating Arteries
    end
    toc % Finish timing artery generation.
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.
    
    disp('%%%%%%% Generating Venous Vessels %%%%%%%%%%%%%%%%%') % Display.
    tic % Start timing vein generation.
    if Option_RRTstar;
        Vessel2 = vesselGenerationRRT_star(Vessel2,Perfusion,GM_WM,VeinGenerationIterations,VeinGenerationWeightFactor,RRTstarEpsilon/VoxelSize); % Generating Veins
    else
        Vessel2 = vesselGenerationRRT(Vessel2,Perfusion,GM_WM,VeinGenerationIterations,VeinGenerationWeightFactor); % Generating Veins
    end
    toc % Finish timing vein generation.
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.
    % Performing all of the vessel generation.
    
    InletPoints = findPoints(InletOriginal,Vessel1);
    OutletPoints = findPoints(OutletOriginal,Vessel2);
    % Inlet and outlet points might have changed values (but not location)
    % so these update those by finding the closest point to the originally
    % defined inlets and outlets.
    
    if Option_Stroke
        StrokeLocation = findPoints(StrokeOriginal,Vessel1);
    end
    % Finding the location of stroke (if used).
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Obtaining Vessel Geometry %%%%%%%%%%%%%%%%
[L1, AL1, Vol1, Davg1, Aavg1] = vesselGeometry(Vessel1,VoxelSize);
[L2, AL2, Vol2, Davg2, Aavg2] = vesselGeometry(Vessel2,VoxelSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Defining Inlet and Outlet Flowrates %%%%%%
InletFlows = InletFlows*BloodFlow;
OutletFlows = OutletFlows*-BloodFlow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%