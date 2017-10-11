% OPTIONAL: Allows overriding of any previously Initialised parameter.
% Script File: overrideInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 20/09/2017
% Description: Allows the option to override any previous values or options
% in the Initialisation files. This is useful if certain set conditions
% should be used repetitively or if a singular conditions should be changed
% over the course of a loop.


%%%%%% Set Overrides  %%%%%%%%%%%%%%%%%%%%%%%%%%%
TrialOverride = 'none';
% Preset overrides for neonate or rat brain. Options can be [none rat
% neonate tracer arm arm2].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Check Overrides  %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(TrialOverride, 'neonate')
    % Models the head scaled down to the size of a neonate.
    
    Option_VoxelSizeAdjustment = true;
    SetVoxelSize = 0.0015 * 11/15; % [m].
    
    Option_TargetPerfusion = true;
    Option_TargetPerfusionConvert = true;
    TargetPerfusion = 30; % [ml/100g/min].
    
    Q_GreyMatter = 0.5*Q_GreyMatter;
    Q_WhiteMatter = 0.5*Q_WhiteMatter;
    
elseif strcmp(TrialOverride, 'rat')
    % Models the brain only scaled down to the size of a rat.
    Option_DomainAdjustment = 'rat';
    
    Option_VoxelSizeAdjustment = true;
    SetVoxelSize = 0.0015 * 0.1; % [m].
    
    Option_LimitedCooling = true;
    LimitedCoolingHeight = 84;
    
    Option_TargetPerfusion = true;
    Option_TargetPerfusionConvert = true;
    TargetPerfusion = 81; % [ml/100g/min].
    
    Option_TemperatureDifference = false;
    T_Out = 23; % [degC].
    H_Out = 43.8; % [W/m2/degC].
    
elseif strcmp(TrialOverride,'tracer')
    % Models a tracer through the brain instead of temperature. The same
    % solving method is used as equations for concentration are analogous
    % to temperature.
    
    Option_DomainAdjustment = 'brainOnly';
    Option_CubeBrain = false;
    
    Kc_b = 1e-15; % Set diffusivity (conductive transfer) to zero.
    Kc_GreyMatter = 1e-15; % Grey Matter thermal conductivity [W/m/degC].
    Kc_WhiteMatter = 1e-15; % White Matter thermal conductivity [W/m/degC].
    
    Cp_b = 1; % Set thermal capacity of blood to one.
    Cp_GreyMatter = 1e-15; % Set thermal capacity of tissue to zero.
    Cp_WhiteMatter = 1e-15; % Set thermal capacity of blood to zero.
    
    Q_GreyMatter = 0; % Grey Matter metabolic rate [W/m3].
    Q_WhiteMatter = 0; % White Matter metabolic rate [W/m3].
    
    Option_GraetzNumber = false;
    
    BloodTemp = 1; % Inlet Concentration set to one.
    
    if Option_CubeBrain
        % For performing tracer test on a linear model.
        CubeSize = 100;
        Option_CubeBrainLinear = true;
        InletTemp = 1;
        Option_TargetPerfusion = true;
        Option_TargetPerfusionConvert = true;
        TargetPerfusion = 100; % [ml/100g/min].
        VesselDiameterToFlowAdjust = 1;
        CubeVoxelSize = 0.003;
    else
        InletTemp = [1 1 1];
    end
    
    H_Out = 0; % No bonductive transfer with environment.
    
    Nu = 0; % No conductive tranfer between vessel and tissue.
    
    Beta23 = Inf; % Infinite transfer with tissue (although no transfer
    % occurs if tissue thermal capacity is set to zero.
    
    Option_TemperatureDifference = false;
    
    Option_TransientSolve = true;
    Timestep = 0.1; % [s]
    TotalTime = 30; % [s]
    SaveTime = 0.1; % [s]
    TransientInitialTemp = 0;
    
    Option_VaryingInlet = true;
    TimeChange = 1; % Time at which change to inlets change [s].
    BloodTempNew = 0; % Temperature of arterial blood used in Pennes Bioheat Equation perfusion term [degC].
    if Option_CubeBrain
        InletTempNew = 0;
    else
        InletTempNew = [0 0 0];
    end
    
elseif strcmp(TrialOverride, 'arm')
    % Attempts to replictate an arm with a cuboid shape.
    Option_CubeBrain = true;
    Option_CubeBrainArm = true;
    CubeSize = 30;
    CubeArmRatio = 8.25;
    CubeVoxelSize = 0.0886/CubeSize;
    
    Arm_NumVessels = 100;
    CubeBrainArmCounterCurrent = false;
    
    Option_BranchTerminationsOnly = true;
    
    Option_AdiabaticBase = false;
    Option_AdiabaticArmEnds = true;
    
    VesselDiameterToFlowAdjust = 1;
    VesselDiameterAdjust = 1;
    
    Option_TargetPerfusion = true;
    Option_TargetPerfusionConvert = false;
    TargetPerfusion = 0.64;
    
    InletTemp = 37;
    
    Option_GraetzNumber = true;
    
    Q_GreyMatter = 320; % Grey Matter metabolic rate [W/m3].
    Q_WhiteMatter = 320; % White Matter metabolic rate [W/m3].
    
    T_Out = 33.5;
    T_OutNew = 10;
    
elseif strcmp(TrialOverride, 'arm2')
    % Creates a linear model for testing vessel/tissue interactions.
    Option_CubeBrain = true;
    Option_CubeBrainArm2 = true;
    CubeSize = 100;
    CubeArmRatio = 1;
    CubeVoxelSize = 0.73/CubeSize;
    
    Arm_NumVessels = 100;
    CubeBrainArmCounterCurrent = false;
    
    Option_GenerateVessels = false;
    Option_LoadPrevious = false;
    
    Option_BranchTerminationsOnly = true;
    
    Option_AdiabaticBase = false;
    Option_AdiabaticArmEnds = true;
    
    VesselDiameterToFlowAdjust = 1;
    VesselDiameterAdjust = 1;
    
    Option_TargetPerfusion = true;
    Option_TargetPerfusionConvert = false;
    TargetPerfusion = 0.64;
    
    InletTemp = 37;
    
    Option_GraetzNumber = false;
    
    Q_GreyMatter = 0; % Grey Matter metabolic rate [W/m3].
    Q_WhiteMatter = 0; % White Matter metabolic rate [W/m3].
    
    T_Out = 36.5;
    T_OutNew = 10;
    
end
