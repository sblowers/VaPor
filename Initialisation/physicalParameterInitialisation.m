% Input Physical Parameters for the simulation.
% Script File: physicalParameterInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Input for all the physical parameters used within the model.

%%%%%% Density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rho_b = 1050; % Blood density [kg/m3].

Rho_GreyMatter = 1030; % Grey Matter density [kg/m3].
Rho_WhiteMatter = 1030; % White Matter density [kg/m3].
Rho_CSFandEyes = 1000; % CSF and Eyes density [kg/m3].
Rho_Skull = 1520; % Skull (bone) density [kg/m3].
Rho_SoftTissue = 1000; % Soft Tissue density [kg/m3].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Specific heat capacity %%%%%%%%%%%%%%%%%%%
Cp_b = 3800; % Blood specificy heat capacity [J/kg/K].

Cp_GreyMatter = 3700; % Grey Matter specific heat capacity [J/kg/degC].
Cp_WhiteMatter = 3700; % White Matter specific heat capacity [J/kg/degC].
Cp_CSFandEyes = 3800; % CSF and Eyes specific heat capacity [J/kg/degC].
Cp_Skull = 2300; % Skull (bone) specific heat capacity [J/kg/degC].
Cp_SoftTissue = 4000; % Soft Tissue specific heat capacity [J/kg/degC].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Thermal Conductivity  %%%%%%%%%%%%%%%%%%%%
Kc_b = 0.5; % Blood thermal conductivity [W/m/degC].

Kc_GreyMatter = 0.49; % Grey Matter thermal conductivity [W/m/degC].
Kc_WhiteMatter = 0.49; % White Matter thermal conductivity [W/m/degC].
Kc_CSFandEyes = 0.5; % CSF and Eyes thermal conductivity [W/m/degC].
Kc_Skull = 1.16; % Skull (bone) thermal conductivity [W/m/degC].
Kc_SoftTissue = 0.342; % Soft Tissue thermal conductivity [W/m/degC].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Metabolism %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_GreyMatter = 16700; % Grey Matter metabolic rate [W/m3].
Q_WhiteMatter = 4175; % White Matter metabolic rate [W/m3].
Q_CSFandEyes = 0; % CSF and Eyes metabolic rate [W/m3].
Q_Skull = 368.3; % Skull (bone) thermal metabolic rate [W/m3].
Q_SoftTissue = 363.4; % Soft Tissue metabolic rate [W/m3].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Predicted Perfusion %%%%%%%%%%%%%%%%%%%%%%
Option_ConvertPerfusion = true;
% If this option is selected (true), then the perfusion values should be
% given in values of [ml/100g/min] which will then be converted into the
% solver units of [kg/m3/s] using the corresponding tissue density. If the
% option is not selected (false) then the perfusion values should be given
% in units [kg/m3/s] directly.

Perfusion_GreyMatter = 80; % Grey Matter predicted perfusion (80ml/100g/min = 14.42kg/m3/s).
Perfusion_WhiteMatter = 20; % White Matter predicted perfusion (20ml/100g/min = 3.605kg/m3/s).
Perfusion_CSFandEyes = 0; % CSF and Eyes predicted perfusion (0ml/100g/min = 0kg/m3/s).
Perfusion_Skull = 1.8; % Skull (bone) thermal predicted perfusion (1.8ml/100g/min = 0.4788kg/m3/s).
Perfusion_SoftTissue = 2.0; % Soft Tissue predicted perfusion (2.0ml/100g/min = 0.35kg/m3/s).
% Units are [ml/100g/min] if Option_ConvertPerfusion = true.
% Units are [kg/m3/s] if Option_ConvertPerfusion = false.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Scale Perfusion %%%%%%%%%%%%%%%%%%%%%%%%%%
Option_TargetPerfusion = false;
% If this option is selected (true), the perfusion values within the brain
% domain (not in surface tissue) will be adjusted to fit the target
% parameter. This option is also affected by 
% Option_TargetPerfusionConvert. If Option_TargetPerfusionConvert is true 
% then units for Target Perfusion must be given in [ml/100g/min] otherwise 
% they should be given in [kg/m3/s].

Option_TargetPerfusionConvert = true;
TargetPerfusion = 30;
% Units are [ml/100g/min] if Option_TargetPerfusionConvert = true. 
% Units are [kg/m3/s] if Option_TargetPerfusionConvert = false.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Porosity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Porosity_GreyMatter = 0.0672; % Grey Matter porosity [].
Porosity_WhiteMatter = 0.0405; % White Matter porosity [].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Viscosity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Visc_b = 0.0035; % Blood viscosity [Pa s].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Capillary Information %%%%%%%%%%%%%%%%%%%%
D_Cap = 10e-6; % Capillary diameter [m].
Tortuosity = 2.5; % Capillary tortuosity [].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



