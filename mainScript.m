% Main Excecution file for the VaPor model.
% Script File: mainScript.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 07/11/2018
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

if ~Option_PennesOnly || (Option_PennesOnly && Option_Stroke)
    % Skip flowrate simulations if only solving Pennes Equation (except
    % when modelling a stroke using Pennes Equation).
    
    vesselReading
    % Creates 1D vessels.
    
    domainIntersection
    % Establishes intersection and mass transfer between 3D domain and 1D
    % vessels.
    
    flowrateSolving1D
    % Solves flowrates in 1D vessels and corresponding diameters.
    
    flowrateSolving3D
    % Solves flowrates in 3D porous domain.
    
    if Option_CounterCurrentFlow
        flowrateSolving3D_CounterCurrent
        % Solves flowrates in 3D porous domain an additional for
        % calculating counter current-flow.
    end
    
    if Option_PseudoCounterCurrentFlow
        U2 = -(Option_PseudoAmount-1)*U; U = Option_PseudoAmount*U; 
        V2 = -(Option_PseudoAmount-1)*V; V = Option_PseudoAmount*V; 
        W2 = -(Option_PseudoAmount-1)*W; W = Option_PseudoAmount*W; 
        % Adjusts the velocities for the use of pseudo counter current
        % flow.
    end
end

if Option_Stroke
    StrokePerfOriginal = MeasuredPerfusion;
    % Old perfusion values measured from velocoties.
    strokeAdjustment
    % Adjusts velocities based on stroke
    StrokePerfNew = MeasuredPerfusion;
    % New perfusion values measured from velocoties.
    
    if Option_PennesOnly
        Perfusion(GM_WM) = Perfusion(GM_WM).*(StrokePerfNew(GM_WM)./StrokePerfOriginal(GM_WM));
        PerfusionAlt(GM_WM) = PerfusionAlt(GM_WM).*(StrokePerfNew(GM_WM)./StrokePerfOriginal(GM_WM));
        % Adjusts perfusion values for use in Pennes Bioheat Equation.
    end
end

if Option_ParticleTracer
    if Option_CounterCurrentFlow
        particleTracerMethod_CounterCurrent
    else
        particleTracerMethod
    end
else
    
    if Option_TransientSolve
        temperatureSolver_CN % Transient solver with Crank-Nicholson timestepping.
    else
        temperatureSolver % Steady-state temperature solver.
        % Solves temperatures across all domains.
        
        if Option_TemperatureDifference
            
            if Option_PennesOnly
                Tt_Save = Tt;
            else
                Tt_Save = Tt;
                Tb_Save = Tb;
                T_Art_Save = T_Art;
                T_Vein_Save = T_Vein;
            end
            % Saving values for comparison after difference is calculated.
            
            H_Out = H_OutNew;
            T_Out = T_OutNew;
            % Set new boundary conditions for second trial.
            
            temperatureSolver
            % Solves temperatures across all domains.
            
            Tt_Diff = Tt - Tt_Save;
            % Establish difference between two results.
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%') % Display.