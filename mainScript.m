% Main Excecution file for the VaPor model.
% Script File: mainScript.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Establishes pathways and excecutes initialisation and 
% excecution scripts in order. 

%%%%%% Cleans Workspace %%%%%%%%%%%%%%%%%%%%%%%%%
clear % Deletes all variables.
close all % Closes all figures.
clc % Clears the command window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%') % Display.

%%%%%% Adds folder pathways for functions %%%%%%%
addPathways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%
physicalParameterInitialisation
% This contains all the input physical parameters for the model.

domainReadingInitialisation
% This contains all the options for reading and creating domains.

vesselReadingInitialisation
% This contains all the options for reading and creating vessels.

temperatureSolverInitialisation
% This contains all the options for solving temperature.

overrideInitialisation
% This contains any overriding options from the previous initialisations.
% This is useful to set up loops that change certain parameters (such as
% increasing number of vessels generated, or different temperatures).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Execution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
domainReading
% Creates 3D domain.

if ~Option_PennesOnly % Skip flowrate simulations if only solving Pennes Equation.
    vesselReading
    % Creates 1D vessels.
    
    domainIntersection
    % Establishes intersection and mass transfer between 3D domain and 1D
    % vessels.
    
    flowrateSolving1D
    % Solves flowrates in 1D vessels and corresponding diameters.
    
    flowrateSolving3D
    % Solves flowrates in 3D porous domain.
end

temperatureSolver
% Solves temperatures across all domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Finding Temperature Difference Option %%%%
if Option_TemperatureDifference 
    
    if Option_PennesOnly
        Tt_Save = Tt;
    else
        Tt_Save = Tt;
        Tb_Save = Tb;
        T_Art_Save = T_Art;
        T_Vein_Save = T_Vein;
    end
    
    H_Out = H_OutNew;
    T_Out = T_OutNew;
    % Set new boundary conditions for second trial.
    
    temperatureSolver
    % Solves temperatures across all domains.
    
    Tt_Diff = Tt - Tt_Save;
    % Establish difference between two results.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%') % Display.