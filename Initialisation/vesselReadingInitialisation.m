% Initialise Information for Implementing 3D Domains
% Script File: vesselReadingInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Initialising parameters for reading and extracting
% the 1D vascular trees. Involves setting file locations for the base
% vascular tree files and what adjustments of the data will be performed
% such as generation of further branches on the tree.

%%%%%% Defining File Locations for Vessels %%%%%%
ArteriesLocation = 'Vessels\';
ArteriesFile = 'BG001.swc';
% File locations for the arteries.

VeinsLocation = 'Vessels\';
VeinsFile = 'Vein_Vessel.mat';
% File locations for the veins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Defining Inlet and Oultets %%%%%%%%%%%%%%%
InletPoints = [1 56 107]; % Points on original vessel tree that are inlets.
OutletPoints = [148 242]; % Points on original vessel tree that are outlets.
InletFlows = [0.2,0.4,0.4]; % Corresponding fraction of overall flowrate as inlet flowrate.
OutletFlows = [0.5,0.5]; % Corresponding fraction of overall flowrate as outlet flowrate.
% The flowrate fractions must sum up to 1 for the vessel flowrates to be
% balanced.

Option_PressureSolve = false;
% Option to solve using pressure however if all inlet and outlet flowrates
% are set then a direct solve is faster and more accurate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Default Adjustments %%%%%%%%%%%
Option_VesselAdjustment = 'default';
% This option adjusts the default vessel trees so that they fit into the
% domain (rotations, scaling and translation etc.). If different trees are
% used, then alterations have to be made to tailor the trees to the domain.
% Currently only ['default'] is defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Loading Previous Vessel Trees %
Option_LoadPrevious = false;
LoadPreviousFilename = 'Previous Results\2500Vessels.mat';
% If a saved results file contains information
% ['Vessel1','Vessel2','InletPoints','OutletPoints'] they will be loaded
% from the given filename instead. This can still be used in conjunction
% with Option_GenerateVessels to expand the vessel tree further by a given
% number of iterations.

Option_LoadPreviousCoarsening = 'none';
% If a previous Vessel tree is to be used but the domain has be coarsened,
% then this option can be used to coarsen the vessel tree accordingly.
% Options can include ['none','two','three','six'].

Option_LoadPreviousRefining = 'none';
% If a previous Vessel tree is to be used but the domain has be refined,
% then this option can be used to refine the vessel tree accordingly.
% Options can include ['none','two','three','six']. (This is effectively
% the opposite of Option_LoadPreviousCoarsening).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Generating Vessels %%%%%%%%%%%%
Option_GenerateVessels = false;
% This option allows generation of vessels inside the brain domain using
% the RRT-Algorithm (Option_GenerateVessels = true).

ArteryGenerationIterations = 20000; % Number of iterations to generate arteries.
VeinGenerationIterations = 20000; % Number of iterations to generate veins.
ArteryGenerationWeightFactor = 2;
VeinGenerationWeightFactor = 2;
% Weight factor alters the weighting function inputted into the RRT solver.
% It adjusts the probablilty by P^WF/sum(P^WF) where P is the local
% probablility and WF is the Weighting Function. It was found that setting
% this to 2 created a better segmentation between grey and white matter but
% beyond that the segmentation becomes too high.

Option_RRTstar = true;
% Option_RRTstar toggles an alternative RRT algorithm whereby all generated
% vessels are limited to a certain length. This gives a more realistic
% branching tree as not all vessels are joined by right-angles.
RRTstarEpsilon = 0.003; % Maximum length of grown vessel [m].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Transfer at Branch Ends %%%%%%%
Option_BranchTerminationsOnly = true;
% If 'true' the transfer of blood between vessel domains and porous domains
% is limited only at the segments that are on branch terminations of vessel
% trees. Otherwise the transfer occurs in all segments. The amount
% transferred is still weighted by lenght of vessel segment in both cases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Simulating Stroke %%%%%%%%%%%%%
Option_Stroke = false;
% Simulates a stroke by obstructing blood down a path on the arterial tree.
% Initially solves the entire flow normally and then reduces upstream and
% downstream flow from the point of obstruction. Overall flow is then NOT
% balanced so it is reduced after this. This option also works if solving
% temperatures using Pennes Bioheat Equation only.

StrokeLocation = 825; % Location of stroke obstruction on arterial tree.
StrokeLocationPos = [825, 3, 41.7308, 46.0255, 20.7829, 0.0034, 824];
% If vessels are loaded then this position is used instead which refers to
% the StrokeLocation = 825 on the original tree, with domain coarsening of
% 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Vessel Diameters %%%%%%%%%%%%%%
VesselDiameterToFlowAdjust = 1; % Adjusts diameters by adjusting flowrate parameter temporarily.
VesselDiameterAdjust = 1; % Adjusts diameters directly.
% Vessel diameters are calculated from a relationship with the flow
% through them. These options can adjust the scale artificially to make
% fatter or thinner vessels overall.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for CubeBrainArm %%%%%%%%%%%%%%%%%%
% If CubeBrainArm is used, these options specify parameters for the vessels
% used.

Arm_NumVessels = 100; % Number of vessel segments along the arm length.
Option_CubeBrainArmCounterCurrent = false;
% If Option_CubeBrainArmCounterCurrent is 'false' a single vessel passes
% through the core of the brain and flow occurs only in one direction. If
% the option is set to 'true' two vessels operating countercurrently along
% the centre of the arm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

