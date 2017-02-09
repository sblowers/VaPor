% Initialise Information for Implementing 3D Domains
% Script File: domainReadingInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
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
% include['none', 'default' or 'rat'].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Adjusting Voxel Size %%%%%%%%%%
Option_VoxelSizeAdjustment = false;
SetVoxelSize = 0.0015; % Desired voxel size [m].
% If a different size voxel is wanted for scaling the model input it
% here and set Option_VoxelSizeAdjustment = true. Default size of voxel is
% 0.0015m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Coarsening Domain Grid %%%%%%%%
Option_DomainCoarsening = 'three';
% The factor for which the grid domain is coarsened in the default
% adjustment option. options include ['none','two','three','six'] as the
% size of the files are divisible by six. The option 'none' is not
% recommended as the size of the domain is too big for the linear solver.
% All published trials were performed with option_coarsening = 'two'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
