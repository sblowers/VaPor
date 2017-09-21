% Initialise Information for Implementing 3D Domains
% Script File: domainReadingInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 20/09/2017
% Description: Initialising parameters for reading and extracting
% the 3D domain. Involves setting file locations for the domain files and
% what adjustments of the data will be performed.


%%%%%% File Locations for Image Sequences %%%%%%%
% File locations for the corresponding image sequences.
GreyMatterLocation = 'Image Sequence\Grey Matter\';
GreyMatterFilePrefix = 'TPM_t001_z';
GreyMatter = zeros(121,145,121); % Initialise size of images and number of images.
GreyMatterDirection = 'z'; % Direction of which images are sliced. ['x' 'y' or 'z'].
% File information for Grey Matter image sequence.

WhiteMatterLocation = 'Image Sequence\White Matter\';
WhiteMatterFilePrefix = 'TPM_t002_z';
WhiteMatter = zeros(121,145,121); % Initialise size of images and number of images.
WhiteMatterDirection = 'z'; % Direction of which images are sliced. ['x' 'y' or 'z'].
% File information for White Matter image sequence.

CSFandEyesLocation = 'Image Sequence\CSF and Eyes\';
CSFandEyesFilePrefix = 'TPM_t003_z';
CSFandEyes = zeros(121,145,121); % Initialise size of images and number of images.
CSFandEyesDirection = 'z'; % Direction of which images are sliced. ['x' 'y' or 'z'].
% File information for CSF and Eyes image sequence.

SkullLocation = 'Image Sequence\Skull\';
SkullFilePrefix = 'TPM_t004_z';
Skull = zeros(121,145,121); % Initialise size of images and number of images.
SkullDirection = 'z'; % Direction of which images are sliced. ['x' 'y' or 'z'].
% File information for Skull image sequence.

SoftTissueLocation = 'Image Sequence\Soft Tissue\';
SoftTissueFilePrefix = 'TPM_t005_z';
SoftTissue = zeros(121,145,121); % Initialise size of images and number of images.
SoftTissueDirection = 'z'; % Direction of which images are sliced. ['x' 'y' or 'z'].
% File information for Soft Tissue image sequence.

MaskLocation = 'Image Sequence\Background Mask\';
MaskFilePrefix = 'TPM_t006_z';
Mask = zeros(121,145,121); % Initialise size of images and number of images.
MaskDirection = 'z'; % Direction of which images are sliced. ['x' 'y' or 'z'].
% File information for Mask image sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Default Adjustments %%%%%%%%%%%
% Options for adjusting the image sequences.
Option_DomainAdjustment = 'default';
% The adjustment parameters are for the default images used. Options
% include['none', 'default', 'brainOnly' or 'rat'].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Adjusting Voxel Size %%%%%%%%%%
Option_VoxelSizeAdjustment = false;
SetVoxelSize = 0.0015; % Desired voxel size [m].
% If a different size voxel is wanted for scaling the model input it
% here and set Option_VoxelSizeAdjustment = true. Default size of voxel is
% 0.0015m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Coarsening Domain Grid %%%%%%%%
Option_DomainCoarsening = 'two';
% The factor for which the grid domain is coarsened in the default
% adjustment option. options include ['none','two','three','six'] as the
% size of the files are divisible by six. The option 'none' is not
% recommended as the size of the domain is too big for the linear solver.
% All published trials were performed with option_coarsening = 'two'.

% Option_DomainCoarsening = true;
% Option_DomainCoarseningFactor = 1/2;
% These options are used for altering the domain by an arbitrary scaling
% value rather than using the above set values. Creates slightly different
% domains compared due to the scaling method.
% CURRENTLY NOT USED: Code is present in domainReading.m but additional
% code is required in vesselReading.m for scaling the vessels accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Counter-Current Porous Flow %%%
Option_CounterCurrentFlow = false;
% Sets whether counter current flow is used in the porous domain. This
% distributes the blood after it leaves the arteries to voxels weighted by
% the perfusion values. A second flow solver is then used to drain the
% blood from these voxels to the veins. Both distributing and draining
% flows are used in the temperature solution.

if Option_CounterCurrentFlow
    Porosity_GreyMatter = Porosity_GreyMatter/2;
    Porosity_WhiteMatter = Porosity_WhiteMatter/2;
    % As there effectively becomes two domains for the counter-current flow
    % in the porous media, they each take up half of the available volume
    % fraction. Therefore the porosity for each volume is halved to account
    % for this.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Pseudo Counter-Current Flow %%%
Option_PseudoCounterCurrentFlow = false;
Option_PseudoAmount = 5; % Pseudo counter-current factor [].
% Instead of modelling a full counter current model, this option just
% splits the net flow between voxels into two parts, forward and reverse.
% The forward flow is multiplied by the factor specified in
% Option_PseudoAmount and the reverse flow is balances the net flow so it
% equals the orginal (it is equal to Forward*-(Option_PseudoAmount-1)).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Cube Brain %%%%%%%%%%%%%%%%%%%%
% These options replace the brain domain with various cubic structures for
% investigation purposes.

Option_CubeBrain = false;
CubeSize = 60; % Number of voxels in cube side [].
CubeVoxelSize = 0.114/CubeSize; % Desired cubic voxel size [m].
% Replaces the domain used with a cubic structure. This is only comprised
% of 50% grey and 50% white matter with no other tissue present. The number
% of voxels in the height width and length of the cube is designated by
% CubeSize and the voxel element size is dictated by CubeVoxel Size.
% This options also specifies that a single artery and vein is created in
% opposite corners of the cube and overrides any other loaded vasculature.
% Generation of vessels can still be applied.
% Option_CubeBrain has to be true for any of the options below to be
% applied.

Option_CubeBrainLinear = false;
% Replaces the domain with a linear structure of voxels (in the
% z-direction) that is has a length of CubeSize voxels but is only one
% voxel wide in both the x and y-direcetions. One artery is connected to
% one end and one vein is connected to the other end.

Option_CubeBrainArm = false;
CubeArmRatio = 8.25; % Ratio of voxels in z-direction to x or y-directions.
% Replaces the domain with a cuboid structure with a vessel passing along
% the centre (approximately representing a forearm). The number of voxels
% in the x and y-directions are given by CubeSize. The number of voxels in
% the z-direction is given by CubeSize*CubeArmRatio (rounded). The tissue
% is comprised of 50% grey and 50% white matter with no other tissue
% present.
% A single or counter-current pair of vessels is modelled (depending on
% options set in vesselReadingInitialisation).

Option_CubeBrainArm2 = false;
% Similar to Option_CubeBrainArm except that the number of voxels in the x
% and y-directions is set to 1. The lenght of the arm stays the same as
% before (CubeSize*CubeArmRatio).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%