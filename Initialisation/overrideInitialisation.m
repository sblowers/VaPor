% OPTIONAL: Allows overriding of any previously Initialised parameter.
% Script File: overrideInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Allows the option to override any previous values or options
% in the Initialisation files. This is useful if certain set conditions
% should be used repetitively or if a singular conditions should be changed
% over the coarse of a loop.


%%%%%% Set Overrides  %%%%%%%%%%%%%%%%%%%%%%%%%%%
TrialOverride = 'none'; 
% Preset overrides for neonate or rat brain. Options can be [none rat
% neonate].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Check Overrides  %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(TrialOverride, 'neonate')
    Option_VoxelSizeAdjustment = true;
    SetVoxelSize = 0.0015 * 11/15; % [m].
    
    Option_TargetPerfusion = true;
    Option_TargetPerfusionConvert = true;
    TargetPerfusion = 30; % [ml/100g/min].
    
    Q_GreyMatter = 0.5*Q_GreyMatter;
    Q_WhiteMatter = 0.5*Q_WhiteMatter;
    
elseif strcmp(TrialOverride, 'rat')
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%