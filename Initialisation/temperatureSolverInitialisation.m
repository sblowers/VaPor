% Input Information for Solving Temperatures
% Script File: temperatureSolverInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Sets the options and boundary conditions for the temperature
% solver.

%%%%%% Setting Blood Boundary Temperatures %%%%%%
BloodTemp = 37; % Temperature of arterial blood used in Pennes Bioheat Equation perfusion term [degC].
InletTemp = [37 37 37]; % Corresponding temperature at inlets [degC].
% Setting the temperature of blood as it enters the model through either
% the inlets of the arterial tree or though the use of Pennes Bioheat
% Equation perfusion term.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Setting Boundary Heat Transfer %%%%%%%%%%%
H_Out = Inf; % Boundary heat transfer coefficient [W/m2/degC].
% This boundary heat transfer coefficient that will be used. If this is set
% to H_Out = Inf then the temperature at all boundary voxels will be set
% to the value for T_Out.

T_Out = 36.5; % Temperature of boundary [degC].
% This is the boundary temperature that will be used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Inter-Domain Heat Transfer in Voxels %%%%%
Beta23 = Inf; % Voxel inter-domain heat transfer coefficient [W/degC].
% The inter-domain heat transfer coefficient for tissue and blood phases
% within the voxels. For finite values, both tissue and blood phases will
% be solved which increases the size of the matrix to be solved. If
% equilibrium between tissue and blood phases is expected, set Beta23 = Inf
% which will solve only a single domain and reduce the memory requirments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Nusselt Number %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nu = 4; % Nusselt Number [].
% Set the Nusselt Number for heat transfer between voxels and vessel
% segments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Option for Pennes Bioheat Equation Only %%
Option_PennesOnly = false;
% If this option is selected (true) then Pennes Bioheat Equation will be
% used to solve the temperatures. Additionally, all parts of the model
% associated with generating flowrates will be skipped (vesselReading,    
% domainIntersection, flowrateSolving1D, flowrateSolving3D).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%% Option for Modelling Temp Difference %%%%%
Option_TemperatureDifference = false;
% This option tells the model to run two trials. The first with the
% boundary conditions [H_Out and T_Out] defined above, and the second with a new
% boundary conditions [H_OutNew and T_OutNew] defined below. This generates
% two sets of results that can be compared quickly.

H_OutNew = Inf; % Boundary heat transfer coefficient [W/m2/degC].
% This is the second boundary heat transfer coefficient that will be used.
% If this is set to H_OutNew = Inf then the temperature at all boundary
% voxels will be set to the value for T_OutNew.

T_OutNew = 10; % Temperature of boundary [degC].
% This is the second boundary temperature that will be used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%% Option for Limited Scalp Cooling %%%%%%%%%
Option_LimitedCooling = false;
% This option sets a height limit for boundary heat transfer. Currently
% only a limit in the Z direction can be implemented through 
% LimitedCoolingHeight. 

LimitedCoolingHeight = 84;
% Any voxels below this voxel value in the Z direction will not interact
% with the boundary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

