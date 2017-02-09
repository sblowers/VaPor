
%%%%%%% Loading Vessels from Files %%%%%%%%%%%%%%
if Option_LoadPrevious
    disp('%%%%%% Loading Previous Results %%%%%%%%%%%%%%%%%%%')
    disp(['Loading: ' LoadPreviousFilename])
    load(LoadPreviousFilename,'Vessel1','Vessel2','InletPoints','OutletPoints')
    disp('Loaded')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    % If Options_LoadPrevious is used then the alterations to the vessel
    % tree do not need to be performed.
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

%%%%%% Coarsening Vessels if Required %%%%%%%%%%%
if Option_LoadPrevious
    if strcmp(Option_LoadPreviousCoarsening, 'two') || strcmp(Option_LoadPreviousCoarsening, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1)./2 + 1; % Scale arterial tree down by a factor of 2
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1)./2 + 1;
    end
    
    if strcmp(Option_LoadPreviousCoarsening, 'three') || strcmp(Option_LoadPreviousCoarsening, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1)./3 + 1; % Scale arterial tree down by a factor of 3
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1)./3 + 1;
    end
    % If Option_LoadPrevious is used then the corresponding coarsening is
    % defined by Option_LoadPreviousCoarsening and not by the domain
    % Option_Coarsening.
    
    if strcmp(Option_LoadPreviousRefining, 'two') || strcmp(Option_LoadPreviousRefining, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1).*2 + 1; % Scale arterial tree down by a factor of 2
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1).*2 + 1;
    end
    
    if strcmp(Option_LoadPreviousRefining, 'three') || strcmp(Option_LoadPreviousRefining, 'six')
        Vessel1(:,3:5) = (Vessel1(:,3:5)-1).*3 + 1; % Scale arterial tree down by a factor of 3
        Vessel2(:,3:5) = (Vessel2(:,3:5)-1).*3 + 1;
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
    % If the domain has been coarsened using Option_Coarsening, then the
    % vessels are scaled to match the new voxel numbers. 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Generating Vessels if Required %%%%%%%%%%%
if Option_GenerateVessels
    InletOriginal = Vessel1(InletPoints,:);
    OutletOriginal = Vessel2(OutletPoints,:);
    % Storing the inlets and outles with their locations so they can be
    % found later.
    
    disp('%%%%%%% Generating Arterial Vessels %%%%%%%%%%%%%%%') % Display.
    tic % Start timing artery generation.
    Vessel1 = vesselGenerationRRT(Vessel1,Perfusion,GM_WM,ArteryGenerationIterations,ArteryGenerationWeightFactor); % Generating Arteries
    toc % Finish timing artery generation.
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.
    
    disp('%%%%%%% Generating Venous Vessels %%%%%%%%%%%%%%%%%') % Display.
    tic % Start timing vein generation.
    Vessel2 = vesselGenerationRRT(Vessel2,Perfusion,GM_WM,VeinGenerationIterations,VeinGenerationWeightFactor); % Generating Veins
    toc % Finish timing vein generation.
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.
    % Performing all of the vessel generation.
    
    InletPoints = findPoints(InletOriginal,Vessel1);
    OutletPoints = findPoints(OutletOriginal,Vessel2);
    % Inlet and outlet points might have changed values (but not location) 
    % so these update those by finding the closest point to the originally
    % defined inlets and outlets.
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
