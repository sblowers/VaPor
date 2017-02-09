% Solves mass flowrates within vessels and calculates vessel diameters.
% Script File: temperatureSolverInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Solves the flowrates in the vessel trees and then calculates
% the diameters of vessels based on the resultant flowrates. As diameter is
% required to solve for flowrates, the process is iterated until a target
% residual value (ResLimit = 1e-10) is reached.


%%%%%% Setting Initial Guess for Diameters %%%%%%
Vessel1(:,6) = 0.001; % Initial guess for Artery Diameters.
Vessel2(:,6) = 0.001; % Initial guess for Vein Diameters.
% Setting an initial guess for the diameters so that they can be solved by
% the solver. This overwrites any diameter information that was present
% when the file is loaded.
[L1, AL1, Vol1, Davg1, Aavg1] = vesselGeometry(Vessel1,VoxelSize);
[L2, AL2, Vol2, Davg2, Aavg2] = vesselGeometry(Vessel2,VoxelSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Solving Flowrates and Diameters %%%%%%%%%%
disp('%%%%%%% Solving 1D Flowrates and Diameters %%%%%%%%') % Display.
tic % Start timing.

ResTest = 0; % Variable used for breaking loop. 
ResLimit = 1e-10; % Set limit for the mean difference in diameters.
iteration = 0; % Counts the number of iterations
while ResTest == 0;
    iteration = iteration + 1; % Increases the iteration number.
    disp(['Iteration: ' num2str(iteration)]); % Displays the current iteration of solving
    
    DartOld = Vessel1(:,6);
    DveinOld = Vessel2(:,6);
    % Saves the diameters for comparison with the derived diameters
    
    [FdotArt, SplitArt] = velocitySolver1DFunction(...
        Vessel1,VoxelSize,MdotArt,Aavg1,L1,Rho_b,Visc_b,InletPoints,InletFlows,[],[]);
    [FdotVein, SplitVein] = velocitySolver1DFunction(...
        Vessel2,VoxelSize,MdotVein,Aavg2,L2,Rho_b,Visc_b,[],[],OutletPoints,OutletFlows);
    % Solves the flowrates within the arteries and veins. Flowrates are
    % generated with two values. The first column refers to
    % the flowrate closer to the connecting node. The second column refers
    % to the flowrate closer to the current node (node N). Positive flows
    % directed towards the current node and negative flows are directed
    % towards the connecting node.
    % Splits in flow occur when the domain mass transfer comes from or  
    % flows towards both sides of the segment. Needs to be checked for 
    % proper transfers in temperature solver. The fraction in the second
    % column refers to the fraction of mass transfer that comes from or
    % flows to the current node (node N).
    
    Dart = 0.0332*abs(FdotArt(:,2)).^0.3703; % Diameter Flowrate correlation.
    Dart(Dart<=1e-6) = 1e-6; % Set minimum diameter as 1micron.
    % Artery diameters based on flowrates. Any diameter less than 1micron
    % is then set to 1micron. This avoids any zero volume line segments at
    % branch terminations.
    
    Dvein = 0.0332*abs(FdotVein(:,2)).^0.3703; % Diameter Flowrate correlation.
    Dvein(Dvein<=1e-6) = 1e-6; % Set minimum diameter as 1micron.
    % Venous diameters based on flowrates. Any diameter less than 1micron
    % is then set to 1micron. This avoids any zero volume line segments at
    % branch terminations.
    
    
    
    Vessel1(:,6) = Dart;
    Vessel2(:,6) = Dvein;
    % Update diameters in vessel tree
    
    if mean(abs(Dart-DartOld)) < ResLimit && mean(abs(Dvein-DveinOld)) < ResLimit
        ResTest = 1;
    end
    % Check to see whether the mean difference in diameters is less than
    % the limit.
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Updating Vessel Geometry %%%%%%%%%%%%%%%%%
[L1, AL1, Vol1, Davg1, Aavg1] = vesselGeometry(Vessel1,VoxelSize);
[L2, AL2, Vol2, Davg2, Aavg2] = vesselGeometry(Vessel2,VoxelSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Fixing Splits %%%%%%%%%%%%%%%%%%%%%%%%%%%%
SplitArt(SplitArt(:,2)==Inf,:) = 0;
SplitVein(SplitVein(:,2)==Inf,:) = 0;
% Removing any splits that have the fraction set as infinity. These can
% occur where flow is zero at one half of the segment but triggers the
% split condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc % Finish timing.
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.
