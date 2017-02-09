% Initialises 3D domain.
% Script File: temperatureSolverInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Extracts the tissue data from image files and alters is
% depending on options defined in domainReadingInitialisation. Physical
% parameters defined in physicalParameterInitialisation are applied based
% on the weighting of tissue values in each voxel. 

%%%%%% Extracting Grey Matter Image Data %%%%%%%%
if     strcmp(GreyMatterDirection, 'x'), M = size(GreyMatter,1);
elseif strcmp(GreyMatterDirection, 'y'), M = size(GreyMatter,2);
elseif strcmp(GreyMatterDirection, 'z'), M = size(GreyMatter,3);
else   error('GreyMatterDirection not ''x'', ''y'', or ''z''')
end
% Determine number of images to scan

for N = 1:M
    if N < 10
        Filename = [GreyMatterLocation GreyMatterFilePrefix '00' num2str(N) '.png']; % File numbers 1-10 begin with '00'
    elseif N < 100
        Filename = [GreyMatterLocation GreyMatterFilePrefix '0' num2str(N) '.png']; % File numbers 11-99 begin with '0'
    else
        Filename = [GreyMatterLocation GreyMatterFilePrefix  num2str(N) '.png']; % File numbers >=100 no prefix
    end
    % Get exact filename for image slice
    
    if     strcmp(GreyMatterDirection, 'x'), GreyMatter(N,:,:)=double(imread(Filename)')/255; 
    elseif strcmp(GreyMatterDirection, 'y'), GreyMatter(:,N,:)=double(imread(Filename)')/255;
    elseif strcmp(GreyMatterDirection, 'z'), GreyMatter(:,:,N)=double(imread(Filename)')/255;
    else   error('GreyMatterDirection not ''x'', ''y'', or ''z''')
    end
    % Read image and convert image data into fraction (between 0 and 1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Extracting White Matter Image Data %%%%%%%
if     strcmp(WhiteMatterDirection, 'x'), M = size(WhiteMatter,1);
elseif strcmp(WhiteMatterDirection, 'y'), M = size(WhiteMatter,2);
elseif strcmp(WhiteMatterDirection, 'z'), M = size(WhiteMatter,3);
else   error('WhiteMatterDirection not ''x'', ''y'', or ''z''')
end
% Determine number of images to scan
    
for N = 1:M
    if N < 10
        Filename = [WhiteMatterLocation WhiteMatterFilePrefix '00' num2str(N) '.png']; % File numbers 1-10 begin with '00'
    elseif N < 100
        Filename = [WhiteMatterLocation WhiteMatterFilePrefix '0' num2str(N) '.png']; % File numbers 11-99 begin with '0'
    else
        Filename = [WhiteMatterLocation WhiteMatterFilePrefix  num2str(N) '.png']; % File numbers >=100 no prefix
    end
    % Get exact filename for image slice
    
    if     strcmp(WhiteMatterDirection, 'x'), WhiteMatter(N,:,:)=double(imread(Filename)')/255; 
    elseif strcmp(WhiteMatterDirection, 'y'), WhiteMatter(:,N,:)=double(imread(Filename)')/255;
    elseif strcmp(WhiteMatterDirection, 'z'), WhiteMatter(:,:,N)=double(imread(Filename)')/255;
    else   error('WhiteMatterDirection not ''x'', ''y'', or ''z''')
    end
    % Read image and convert image data into fraction (between 0 and 1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Extracting CSF and Eyes Image Data %%%%%%%
if     strcmp(CSFandEyesDirection, 'x'), M = size(CSFandEyes,1);
elseif strcmp(CSFandEyesDirection, 'y'), M = size(CSFandEyes,2);
elseif strcmp(CSFandEyesDirection, 'z'), M = size(CSFandEyes,3);
else   error('CSFandEyesMatterDirection not ''x'', ''y'', or ''z''')
end
% Determine number of images to scan
    
for N = 1:M
    if N < 10
        Filename = [CSFandEyesLocation CSFandEyesFilePrefix '00' num2str(N) '.png']; % File numbers 1-10 begin with '00'
    elseif N < 100
        Filename = [CSFandEyesLocation CSFandEyesFilePrefix '0' num2str(N) '.png']; % File numbers 11-99 begin with '0'
    else
        Filename = [CSFandEyesLocation CSFandEyesFilePrefix  num2str(N) '.png']; % File numbers >=100 no prefix
    end
    % Get exact filename for image slice
    
    if     strcmp(CSFandEyesDirection, 'x'), CSFandEyes(N,:,:)=double(imread(Filename)')/255; 
    elseif strcmp(CSFandEyesDirection, 'y'), CSFandEyes(:,N,:)=double(imread(Filename)')/255;
    elseif strcmp(CSFandEyesDirection, 'z'), CSFandEyes(:,:,N)=double(imread(Filename)')/255;
    else   error('CSFandEyesDirection not ''x'', ''y'', or ''z''')
    end
    % Read image and convert image data into fraction (between 0 and 1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Extracting Skull Image Data %%%%%%%%%%%%%%
if     strcmp(SkullDirection, 'x'), M = size(Skull,1);
elseif strcmp(SkullDirection, 'y'), M = size(Skull,2);
elseif strcmp(SkullDirection, 'z'), M = size(Skull,3);
else   error('SkullDirection not ''x'', ''y'', or ''z''')
end
% Determine number of images to scan
    
for N = 1:M
    if N < 10
        Filename = [SkullLocation SkullFilePrefix '00' num2str(N) '.png']; % File numbers 1-10 begin with '00'
    elseif N < 100
        Filename = [SkullLocation SkullFilePrefix '0' num2str(N) '.png']; % File numbers 11-99 begin with '0'
    else
        Filename = [SkullLocation SkullFilePrefix  num2str(N) '.png']; % File numbers >=100 no prefix
    end
    % Get exact filename for image slice
    
    if     strcmp(SkullDirection, 'x'), Skull(N,:,:)=double(imread(Filename)')/255; 
    elseif strcmp(SkullDirection, 'y'), Skull(:,N,:)=double(imread(Filename)')/255;
    elseif strcmp(SkullDirection, 'z'), Skull(:,:,N)=double(imread(Filename)')/255;
    else   error('SkullDirection not ''x'', ''y'', or ''z''')
    end
    % Read image and convert image data into fraction (between 0 and 1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Extracting Soft Tissue Image Data %%%%%%%%
if     strcmp(SoftTissueDirection, 'x'), M = size(SoftTissue,1);
elseif strcmp(SoftTissueDirection, 'y'), M = size(SoftTissue,2);
elseif strcmp(SoftTissueDirection, 'z'), M = size(SoftTissue,3);
else   error('SoftTissueDirection not ''x'', ''y'', or ''z''')
end
% Determine number of images to scan
    
for N = 1:M
    if N < 10
        Filename = [SoftTissueLocation SoftTissueFilePrefix '00' num2str(N) '.png']; % File numbers 1-10 begin with '00'
    elseif N < 100
        Filename = [SoftTissueLocation SoftTissueFilePrefix '0' num2str(N) '.png']; % File numbers 11-99 begin with '0'
    else
        Filename = [SoftTissueLocation SoftTissueFilePrefix  num2str(N) '.png']; % File numbers >=100 no prefix
    end
    % Get exact filename for image slice
    
    if     strcmp(SoftTissueDirection, 'x'), SoftTissue(N,:,:)=double(imread(Filename)')/255; 
    elseif strcmp(SoftTissueDirection, 'y'), SoftTissue(:,N,:)=double(imread(Filename)')/255;
    elseif strcmp(SoftTissueDirection, 'z'), SoftTissue(:,:,N)=double(imread(Filename)')/255;
    else   error('SoftTissueDirection not ''x'', ''y'', or ''z''')
    end
    % Read image and convert image data into fraction (between 0 and 1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Extracting Background Mask Image Data %%%%
if     strcmp(MaskDirection, 'x'), M = size(Mask,1);
elseif strcmp(MaskDirection, 'y'), M = size(Mask,2);
elseif strcmp(MaskDirection, 'z'), M = size(Mask,3);
else   error('Mask not ''x'', ''y'', or ''z''')
end
% Determine number of images to scan
    
for N = 1:M
    if N < 10
        Filename = [MaskLocation MaskFilePrefix '00' num2str(N) '.png']; % File numbers 1-10 begin with '00'
    elseif N < 100
        Filename = [MaskLocation MaskFilePrefix '0' num2str(N) '.png']; % File numbers 11-99 begin with '0'
    else
        Filename = [MaskLocation MaskFilePrefix  num2str(N) '.png']; % File numbers >=100 no prefix
    end
    % Get exact filename for image slice
    
    if     strcmp(MaskDirection, 'x'), Mask(N,:,:)=double(imread(Filename)')/255; 
    elseif strcmp(MaskDirection, 'y'), Mask(:,N,:)=double(imread(Filename)')/255;
    elseif strcmp(MaskDirection, 'z'), Mask(:,:,N)=double(imread(Filename)')/255;
    else   error('MaskDirection not ''x'', ''y'', or ''z''')
    end
    % Read image and convert image data into fraction (between 0 and 1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Default Option Adjustments %%%%%%%%%%%%%%%
if strcmp(Option_DomainAdjustment, 'default') || strcmp(Option_DomainAdjustment, 'rat')
    GreyMatter = GreyMatter(1:end-1,1:end-1,1:end-1);
    WhiteMatter = WhiteMatter(1:end-1,1:end-1,1:end-1);
    CSFandEyes = CSFandEyes(1:end-1,1:end-1,1:end-1);
    Skull = Skull(1:end-1,1:end-1,1:end-1);
    SoftTissue = SoftTissue(1:end-1,1:end-1,1:end-1);
    Mask = Mask(1:end-1,1:end-1,1:end-1);
    % Reduces the size of domain by 1 in all directions so that they are
    % all divisible by 6. These are empty slices so no data is lost.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Setting Voxel Size %%%%%%%%%%%%%%%%%%%%%%%
if Option_VoxelSizeAdjustment
    VoxelSize = SetVoxelSize;
else
    VoxelSize = 0.0015;
end
% Set Voxel size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Coarsening Grids %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Option_DomainCoarsening, 'two') || strcmp(Option_DomainCoarsening, 'six')
    GreyMatter = coarsenCubeSize2(GreyMatter); % Coarsens the grid using 2x2x2 cube
    WhiteMatter = coarsenCubeSize2(WhiteMatter); % Coarsens the grid using 2x2x2 cube
    CSFandEyes = coarsenCubeSize2(CSFandEyes); % Coarsens the grid using 2x2x2 cube
    Skull = coarsenCubeSize2(Skull); % Coarsens the grid using 2x2x2 cube
    SoftTissue = coarsenCubeSize2(SoftTissue); % Coarsens the grid using 2x2x2 cube
    Mask = coarsenCubeSize2(Mask); % Coarsens the grid using 2x2x2 cube
    VoxelSize = 2*VoxelSize; % from 2x2x2 coarsening
end
% Coarsening the grid by taking an average of values within a 2x2x2 kernal.
% This halves the size of the domain.

if strcmp(Option_DomainCoarsening, 'three') || strcmp(Option_DomainCoarsening, 'six')
    GreyMatter = coarsenCubeSize3(GreyMatter); % Coarsens the grid using 3x3x3 cube
    WhiteMatter = coarsenCubeSize3(WhiteMatter); % Coarsens the grid using 3x3x3 cube
    CSFandEyes = coarsenCubeSize3(CSFandEyes); % Coarsens the grid using 3x3x3 cube
    Skull = coarsenCubeSize3(Skull); % Coarsens the grid using 3x3x3 cube
    SoftTissue = coarsenCubeSize3(SoftTissue); % Coarsens the grid using 3x3x3 cube
    Mask = coarsenCubeSize3(Mask); % Coarsens the grid using 3x3x3 cube
    VoxelSize = 3*VoxelSize; % from 3x3x3 coarsening
end
% Coarsening the grid by taking an average of values within a 3x3x3 kernal.
% This reduces the size of the domain by a third.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Miscellaneous Adjustments of Domains %%%%%
if strcmp(Option_DomainAdjustment, 'default') || strcmp(Option_DomainAdjustment, 'rat')
    % Fixing Soft Tissue Stuff
    SoftTissue(SoftTissue < 0.04) = 0;
    % Reduces the amount of soft tissue within the domain. This makes the
    % domain size more realistic as it culls some fringe voxels.
    
    % Filling in blank voxels
    boolean = ~logical(Mask) & ~logical(GreyMatter) & ~logical(WhiteMatter) & ~logical(CSFandEyes) & ~logical(Skull) & ~logical(SoftTissue);
    SoftTissue(boolean) = 1; 
    % If any voxel for any reason is blank, it is filled with Soft Tissue
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Defining Domains %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Option_DomainAdjustment, 'default') || strcmp(Option_DomainAdjustment, 'rat')
    % Reducing the size of brain domain
    Cull = 1; %fraction of gives GM_WM volume of around 1500 cm^3
    GM_WM = GreyMatter + WhiteMatter;
    GM_WM(GM_WM<Cull) = 0; % Removes voxels with less than cull amount.
    GM_WM = logical(GM_WM); % Create identifier for brain domain so it can be easily referenced to.
    % Culling certain amount of Grey and White matter from fringe voxels to
    % achieve a more reasonable volume of brain.
    
    % Ensuring only a single brain domain exists
    CC = bwconncomp(GM_WM,6); % Find connectivity of GM_WM domain
    NumPixels = cellfun(@numel,CC.PixelIdxList); % Count number of voxels
    [Biggest,Index] = max(NumPixels); % Find biggest connected volume
    GM_WM = false(size(GM_WM));
    GM_WM(CC.PixelIdxList{Index}) = true; % Only voxels in largest domain remain in GM_WM
    CC = bwconncomp(GM_WM,6); % CC.NumbObjects should now equal 1.
    % If the brain domain is split up at tall, this removes the smaller
    % parts. The removed parts are then converted into Soft Tissue.
    
    % Converting tissue outside of brain region into soft tissue.
    SoftTissue(~GM_WM) = SoftTissue(~GM_WM) + GreyMatter(~GM_WM) + WhiteMatter(~GM_WM);
    GreyMatter(~GM_WM) = 0;
    WhiteMatter(~GM_WM) = 0;
    % Any grey or white matter tissue outside GM_WM domain is converted to
    % soft tissue.
    
    CSFandEyes(GM_WM) = 0;
    Skull(GM_WM) = 0;
    SoftTissue(GM_WM) = 0;
    % Any other tissue in brain domain is removed
    
else
    GM_WM = logical(GreyMatter + WhiteMatter);
    % Creates identifier for brain domain so it can be easily referenced in
    % future.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Overwriting if Rat Brain is Used %%%%%%%%%
if strcmp(Option_DomainAdjustment,  'rat');
    % Rat Brain Options
    CSFandEyes = zeros(size(CSFandEyes));
    Skull = zeros(size(Skull));
    SoftTissue = zeros(size(SoftTissue));
    % Only the brain remains. Surrounding tissue is removed.

    fillGM_WM % Fills the gap inside domain left by removing tissue.
    GreyMatter(GM_WM_filled & ~GM_WM) = 0;
    WhiteMatter(GM_WM_filled & ~GM_WM) = 1;
    % All filled tissue is considered White Matter
    GM_WM = GM_WM_filled;% Updates domain to include filled region.
    % This fills any space left from removing tissue to white matter so there
    % no empty spaces when solving.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Defining Other Domains Used %%%%%%%%%%%%%%
GM = GreyMatter./(GreyMatter + WhiteMatter);
GM = logical(GM >= 0.5);
WM = logical(GM_WM & ~GM);
% Creating identifiers for grey and white matter independently so they can
% be easily referenced in future.

DomTot = logical(GreyMatter + WhiteMatter + CSFandEyes + Skull + SoftTissue);
% Creates identifier for whole domain so it can be easily referenced in
% future.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Converting Perfusion if Required %%%%%%%%%
if Option_ConvertPerfusion 
    Perfusion_GreyMatter = Perfusion_GreyMatter/60*10*Rho_GreyMatter*Rho_b/10^6; % kg/m3/s; % Grey Matter predicted perfusion [kg/m3/s]
    Perfusion_WhiteMatter = Perfusion_WhiteMatter/60*10*Rho_WhiteMatter*Rho_b/10^6; % White Matter predicted perfusion [kg/m3/s]
    Perfusion_CSFandEyes = Perfusion_CSFandEyes/60*10*Rho_CSFandEyes*Rho_b/10^6; % CSF and Eyes predicted perfusion [kg/m3/s]
    Perfusion_Skull = Perfusion_Skull/60*10*Rho_Skull*Rho_b/10^6; % Skull (bone) thermal predicted perfusion [kg/m3/s]
    Perfusion_SoftTissue = Perfusion_SoftTissue/60*10*Rho_SoftTissue*Rho_b/10^6; % Soft Tissue predicted perfusion [kg/m3/s]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Defining Domain specific Parameters %%%%%%
Rho = zeros(size(DomTot));
Rho(DomTot) = (GreyMatter(DomTot)*Rho_GreyMatter + WhiteMatter(DomTot)*Rho_WhiteMatter...
    + CSFandEyes(DomTot)*Rho_CSFandEyes + Skull(DomTot)*Rho_Skull + SoftTissue(DomTot)*Rho_SoftTissue)...
    ./(GreyMatter(DomTot) + WhiteMatter(DomTot) + CSFandEyes(DomTot) + Skull(DomTot) + SoftTissue(DomTot));
% Define voxel tissue density.

Cp = zeros(size(DomTot));
Cp(DomTot) = (GreyMatter(DomTot)*Cp_GreyMatter + WhiteMatter(DomTot)*Cp_WhiteMatter...
    + CSFandEyes(DomTot)*Cp_CSFandEyes + Skull(DomTot)*Cp_Skull + SoftTissue(DomTot)*Cp_SoftTissue)...
    ./(GreyMatter(DomTot) + WhiteMatter(DomTot) + CSFandEyes(DomTot) + Skull(DomTot) + SoftTissue(DomTot));
% Define voxel tissue specific heat capacity.

Kc = zeros(size(DomTot));
Kc(DomTot) = (GreyMatter(DomTot)*Kc_GreyMatter + WhiteMatter(DomTot)*Kc_WhiteMatter...
    + CSFandEyes(DomTot)*Kc_CSFandEyes + Skull(DomTot)*Kc_Skull + SoftTissue(DomTot)*Kc_SoftTissue)...
    ./(GreyMatter(DomTot) + WhiteMatter(DomTot) + CSFandEyes(DomTot) + Skull(DomTot) + SoftTissue(DomTot));
% Define voxel tissue thermal conductivity.

Q = zeros(size(DomTot));
Q(DomTot) = (GreyMatter(DomTot)*Q_GreyMatter + WhiteMatter(DomTot)*Q_WhiteMatter...
    + CSFandEyes(DomTot)*Q_CSFandEyes + Skull(DomTot)*Q_Skull + SoftTissue(DomTot)*Q_SoftTissue)...
    ./(GreyMatter(DomTot) + WhiteMatter(DomTot) + CSFandEyes(DomTot) + Skull(DomTot) + SoftTissue(DomTot));
% Define voxel tissue metabolic heat rate.

Perfusion = zeros(size(DomTot));
Perfusion(DomTot) = (GreyMatter(DomTot)*Perfusion_GreyMatter + WhiteMatter(DomTot)*Perfusion_WhiteMatter...
    + CSFandEyes(DomTot)*Perfusion_CSFandEyes + Skull(DomTot)*Perfusion_Skull + SoftTissue(DomTot)*Perfusion_SoftTissue)...
    ./(GreyMatter(DomTot) + WhiteMatter(DomTot) + CSFandEyes(DomTot) + Skull(DomTot) + SoftTissue(DomTot));
% Define voxel tissue perfusion rate. Used in both temperature solving for
% tissue outside the brain domain and in creating the vessel tree.

Porosity = zeros(size(DomTot));
Porosity(GM_WM) = (GreyMatter(GM_WM)*Porosity_GreyMatter + WhiteMatter(GM_WM)*Porosity_WhiteMatter)./...
(GreyMatter(GM_WM) + WhiteMatter(GM_WM));
% Define voxel tissue porosity. Voxels outside of brain domain have zero
% porosity.

Rho(~DomTot) = NaN;
Cp(~DomTot) = NaN;
Kc(~DomTot) = NaN;
Q(~DomTot) = NaN;
Perfusion(~DomTot) = NaN;
Porosity(~DomTot) = NaN;
% Setting all values outside of the domain to NaN. This facilitates
% visualisation of data with 'planecut'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Establishing Derived Variables %%%%%%%%%%%
DomainVolume = sum(GM_WM(GM_WM))*VoxelSize^3; % Calculating brain volume [m3].

Perf = mean(Perfusion(GM_WM)); % Caluclating mean brain perfusion value [kg/m3/s].
PerfusionAlt = Perfusion;
PerfusionAlt(DomTot) = PerfusionAlt(DomTot)*60/10./Rho(DomTot)/Rho_b*10^6; % Alternative units for Perfusion [ml/100g/min].
PerfAlt = mean(PerfusionAlt(GM_WM)); % Caluclating mean brain perfusion value in alternative units [ml/100g/min].
BloodFlow = DomainVolume*Perf; % Overall brain inlet blood flow [kg/s].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Scaling Perfusion if Required %%%%%%%%%%%%
if Option_TargetPerfusion
    if Option_TargetPerfusionConvert 
        Perfusion(GM_WM) = Perfusion(GM_WM)*(TargetPerfusion/mean(PerfusionAlt(GM_WM)));
    else
        Perfusion(GM_WM) = Perfusion(GM_WM)*(TargetPerfusion/mean(Perfusion(GM_WM)));
    end

    Perf = mean(Perfusion(GM_WM)); % Caluclating mean brain perfusion value [kg/m3/s].
    PerfusionAlt = Perfusion;
    PerfusionAlt(DomTot) = PerfusionAlt(DomTot)*60/10./Rho(DomTot)/Rho_b*10^6; % Alternative units for Perfusion [ml/100g/min].
    PerfAlt = mean(PerfusionAlt(GM_WM)); % Caluclating mean brain perfusion value in alternative units [ml/100g/min].
    BloodFlow = DomainVolume*Perf; % Overall brain inlet blood flow [kg/s].
    % Updating Derived Perfusion Variables

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Finding Domain Borders %%%%%%%%%%%%%%%%%%%
Borders = bwperim(DomTot);
% Establishing node data and boundaries for Temperature Solver.

BordersScalp = zeros(size(DomTot)+[0,0,1]);
BordersScalp(:,:,1) = DomTot(:,:,1);
BordersScalp(:,:,2:end) = DomTot;
BordersScalp = bwperim(BordersScalp);
BordersScalp = BordersScalp(:,:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



