% Solving the flowrates in the 3D porous region.


disp('%%%%%%% Solving 3D Flowrates %%%%%%%%%%%%%%%%%%%%%%') % Display.


%%%%%% Setting up Reference Lists for GM_WM %%%%%
row_convert = find(GM_WM);  % goes from row number to I,J,K
GM_WM_convert = zeros(size(GM_WM));
GM_WM_convert(row_convert) = 1:numel(GM_WM(GM_WM)); % goes from I,J,K to row number
% The sparse matrix to solve only needs those voxels that are within GM_WM.
% These lists help convert sparse matrix row number to location in GM_WM
% (given by I,J,K) and from locations in GM_WM to row numbers in the sparse
% matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initialising Storage Matrices %%%%%%%%%%%%
P_DataTemp = [];
P_DataRowcount = 1;
P_Super = zeros(1e8,3);
D = zeros(numel(GM_WM(GM_WM)),1);
% P_super will be a list of all interactions that will be fed into the
% sparse matrix stored as [row, col, value]. To create the list, it is fed
% by P_DataTemp on every iteration through the domain which is filled up
% by interactions. P_DataRowcount keeps track of where to add P_DataTemp
% at the end of P_Super.
% D contains all the values that are not variable dependent (in this case
% all the Mdot values). This will be the RHS of the equation in the linear
% solver (Mx = D where M is the sparse matrix, x is the list of variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Interactions %%%%%%%%%%%%%%%%%%%
disp('Compiling Interactions within Flow Domain') % Display.
tic % Start timing.

Check = 0; % Used as a condition to set the first voxel to P=0.
for I = 1:size(GM_WM,1)
    for J = 1:size(GM_WM,2)
        for K = 1:size(GM_WM,3)
            if GM_WM(I,J,K)
                % For all voxels, if they are within GM_WM.
                
                %%%%%% Set Pressure Condition %%%%%%%%%%%%%%%%%%%
                if Check == 0
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),1]; % Set P(I,J,K)=0;
                    Check = 1; % Changes condition.
                end
                % Sets the first value to P = 0 so that the solution is
                % bounded.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Mass Source Term %%%%%%%%%%%%%%%%%%%%%%%%%
                D(GM_WM_convert(I,J,K)) = D(GM_WM_convert(I,J,K)) + MdotB(I,J,K);
                % Source term from inter-domain mass transfer.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Interactions with Neighbouring Voxels %%%%
                if I > 1 && GM_WM(I-1,J,K)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I-1,J,K)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I-1,J,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at I-1.
                
                if I < size(GM_WM,1) && GM_WM(I+1,J,K)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I+1,J,K)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I+1,J,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at I+1.
                
                if J > 1 && GM_WM(I,J-1,K)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I,J-1,K)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J-1,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at J-1.
                
                if J < size(GM_WM,2) && GM_WM(I,J+1,K)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I,J+1,K)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J+1,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at J+1.
                
                if K > 1 && GM_WM(I,J,K-1)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I,J,K-1)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K-1),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at K-1.
                
                if K < size(GM_WM,3) && GM_WM(I,J,K+1)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I,J,K+1)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K+1),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;GM_WM_convert(I,J,K),GM_WM_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at K+1.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
                P_Super(P_DataRowcount:P_DataRowcount-1+size(P_DataTemp,1),:) = P_DataTemp; % Add P_DataTemp to P_Super.
                P_DataRowcount = P_DataRowcount + size(P_DataTemp,1); % Update P_DataRowcount.
                P_DataTemp = []; % Clear P_DataTemp.
                % Compiling interactions into P_Super list.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
            end
        end
    end
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Sparse Matrix %%%%%%%%%%%%%%%%%%%
P_Boolean = logical(P_Super(:,1));
P_Super = P_Super(P_Boolean,:); % Delete any excess rows from P_Super.
P_Solve = sparse(P_Super(:,1),P_Super(:,2),P_Super(:,3)); % Convert P_Super into sparse matrix.
clearvars P_super P_Boolean % Delete variables to free up memory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Solving Linear System %%%%%%%%%%%%%%%%%%%%
disp('Starting Linear Solver') % Display.
disp(['Solving Matrix Inversion size: ' num2str(size(P_Solve,1)) ' by ' num2str(size(P_Solve,2))]) % Display.
tic % Start timing.

P_N = P_Solve\D; % Solving the linear system

toc % Finish timing.
disp('Matrix Inversion Completed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Reshaping Resutls from Linear Solver %%%%%
P = zeros(size(GM_WM));
P(row_convert) = P_N;
P(~GM_WM) = NaN;
P = P-min(min(min(P)));
% Converting the results back into the same shape as GM_WM domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Deriving Velocities %%%%%%%%%%%%%%%%%%%%%%
U2=zeros(size(GM_WM,1)+1,size(GM_WM,2),size(GM_WM,3));
V2=zeros(size(GM_WM,1),size(GM_WM,2)+1,size(GM_WM,3));
W2=zeros(size(GM_WM,1),size(GM_WM,2),size(GM_WM,3)+1);
% Initialise velocity matrices.

for I = 1:size(GM_WM,1)
    for J = 1:size(GM_WM,2)
        for K = 1:size(GM_WM,3)
            if GM_WM(I,J,K)
                % For all voxels, if they are within GM_WM.
                
                if I < size(GM_WM,1) && GM_WM(I+1,J,K)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I+1,J,K)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    U2(I+1,J,K) = 1/Rho_b*1/VoxelSize^2*G*(P(I,J,K)-P(I+1,J,K)); % Calculate velocity.
                end
                % create U velocity at voxel boundary I-1.
                
                if J < size(GM_WM,2) && GM_WM(I,J+1,K)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I,J+1,K)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    V2(I,J+1,K) = 1/Rho_b*1/VoxelSize^2*G*(P(I,J,K)-P(I,J+1,K)); % Calculate velocity.
                end
                % create V velocity at voxel boundary J-1.
                
                if K < size(GM_WM,3) && GM_WM(I,J,K+1)
                    AvgPorosity = 0.5*(Porosity(I,J,K) + Porosity(I,J,K+1)); % Average porosity of two voxels.
                    G = Rho_b*VoxelSize*AvgPorosity*pi*D_Cap^2/(32*Visc_b*Tortuosity); % Calculate conductance.
                    W2(I,J,K+1) = 1/Rho_b*1/VoxelSize^2*G*(P(I,J,K)-P(I,J,K+1)); % Calculate velocity.
                end
                % create W velocity at voxel boundary K-1.
                
            end
        end
    end
end
% Velocities are derived from the pressure differences of each voxel. Only
% the faces at I-1, J-1 & K-1 need to be done for each voxel. Velocities at
% boudries are zero by default.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Calculate Derived Perfusion %%%%%%%%%%%%%%
measurePerfusion
% This calculates perfusion values from flows into every voxel. Allows
% comparison with predicted Perfusion values inputted into the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Colocated Velocities %%%%%%%%%%%%
UU2 = 0.5*(U2(1:size(GM_WM,1),1:size(GM_WM,2),1:size(GM_WM,3))+U2(2:size(GM_WM,1)+1,1:size(GM_WM,2),1:size(GM_WM,3)));
VV2 = 0.5*(V2(1:size(GM_WM,1),1:size(GM_WM,2),1:size(GM_WM,3))+V2(1:size(GM_WM,1),2:size(GM_WM,2)+1,1:size(GM_WM,3)));
WW2 = 0.5*(W2(1:size(GM_WM,1),1:size(GM_WM,2),1:size(GM_WM,3))+W2(1:size(GM_WM,1),1:size(GM_WM,2),2:size(GM_WM,3)+1));
Velocity2 = sqrt(UU2.^2 + VV2.^2 + WW2.^2);
% The velocities used in the solver are stored at voxel boundaries which
% are difficult to visualise alongside data stored at voxel centres.
% Therefore combining the average values of two faces gives an
% approximation for values at the centre. These are only used for display
% purposes.

UU2(~GM_WM) = NaN;
VV2(~GM_WM) = NaN;
WW2(~GM_WM) = NaN;
Velocity2(~GM_WM) = NaN;
% Setting all values outside of the domain to NaN. This facilitates
% visualisation of data with 'planecut'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.
