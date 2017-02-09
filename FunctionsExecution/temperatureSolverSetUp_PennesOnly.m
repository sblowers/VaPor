
%%%%%% Setting Up Reference Lists for DomTot %%%%
NumDomTot = numel(DomTot(DomTot));
RowConvert = find(DomTot);  % goes from row number to I,J,K
DomTotConvert = zeros(size(DomTot));
DomTotConvert(RowConvert) = 1:NumDomTot;  % goes from I,J,K to row number
% The sparse matrix to solve only needs those voxels that are within
% DomTot. These lists help convert sparse matrix row number to location in
% DomTot (given by I,J,K) and from locations in DomTot to row numbers in
% the sparse matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initialising Storage Matrices %%%%%%%%%%%%
T_DataTemp = [];
T_DataRowcount = 1;
T_Super = zeros(1e8,3);
D = zeros(NumDomTot,1);
% T_super will be a list of all interactions that will be fed into the
% sparse matrix stored as [row, col, value]. To create the list, it is fed
% by T_DataTemp on every iteration through the domain which is filled up
% by interactions. T_DataRowcount keeps track of where to add T_DataTemp
% at the end of T_Super.
% D contains all the values that are not variable dependent (in this case
% all the Mdot values). This will be the RHS of the equation in the linear
% solver (Mx = D where M is the sparse matrix, x is the list of variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling 3D Voxel Interactions %%%%%%%%%%
disp('Compiling Interactions within 3D Voxels') % Display.
tic % Start timing.

for I=1:size(DomTot,1)
    for J=1:size(DomTot,2)
        for K=1:size(DomTot,3)
            if DomTot(I,J,K)
                % For all voxels, if they are within DomTot.
                
                
                %%%%%% Establishing Row Values %%%%%%%%%%%%%%%%%%
                Row = DomTotConvert(I,J,K);
                if I<size(GM_WM,1), RowUI = DomTotConvert(I+1,J,K); end
                if I>1,             RowDI = DomTotConvert(I-1,J,K); end
                if J<size(GM_WM,2), RowUJ = DomTotConvert(I,J+1,K); end
                if J>1,             RowDJ = DomTotConvert(I,J-1,K); end
                if K<size(GM_WM,3), RowUK = DomTotConvert(I,J,K+1); end
                if K>1,             RowDK = DomTotConvert(I,J,K-1); end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Preallocating Searchable Parameters %%%%%%
                Vol=VoxelSize^3;
                Ax=VoxelSize^2;
                Ay=VoxelSize^2;
                Az=VoxelSize^2;
                % Volume and areas of voxels
                
                Rho_t = Rho(I,J,K);
                Cp_t = Cp(I,J,K);
                Kc_t = Kc(I,J,K);
                Q_t = Q(I,J,K);
                % Establishing searchable parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Boundary Interactions %%%%%%%%%%%%%%%%%%%%
                if Borders(I,J,K) % If voxel lies on domain boundary.
                    BorderMinX=0; BorderMaxX=0; BorderMinY=0; BorderMaxY=0; BorderMinZ=0; BorderMaxZ=0;
                    % Sets all boundary heat transfer to zero
                    
                    
                    %%%%%% Check for Boundary Heat Transfer %%%%%
                    if ~Option_LimitedCooling || K >= LimitedCoolingHeight % Check if Option_LimitedCooling is disabled or both enabled and conditions are met.
                        if I==1              || ~DomTot(I-1,J,K), BorderMinX = 1; end
                        if I==size(DomTot,1) || ~DomTot(I+1,J,K), BorderMaxX = 1; end
                        if J==1              || ~DomTot(I,J-1,K), BorderMinY = 1; end
                        if J==size(DomTot,2) || ~DomTot(I,J+1,K), BorderMaxY = 1; end
                        if K==1              || ~DomTot(I,J,K-1), BorderMinZ = 1; end
                        if K==size(DomTot,3) || ~DomTot(I,J,K+1), BorderMaxZ = 1; end
                        % Check which faces are on bounderies for heat
                        % transfer to occur.
                    end
                    
                    if K == 1, BorderMinZ=0; end
                    % Set base of model to adiabatic.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Boundary Interactions %%%%%%%%%%
                    if H_Out == Inf % If boundary heat transfer is infinite, set voxel temperature equal to boundary temperature.
                        D(Row) = D(Row) + -(BorderMinX+BorderMaxX)*T_Out; % Tissue domain boundaries in the X direction.
                        D(Row) = D(Row) + -(BorderMinY+BorderMaxY)*T_Out; % Tissue domain boundaries in the Y direction.
                        D(Row) = D(Row) + -(BorderMinZ+BorderMaxZ)*T_Out; % Tissue domain boundaries in the Z direction.
                        T_DataTemp = [T_DataTemp;Row,Row,-(BorderMinX+BorderMaxX)]; % Tissue domain boundaries in the X direction.
                        T_DataTemp = [T_DataTemp;Row,Row,-(BorderMinY+BorderMaxY)]; % Tissue domain boundaries in the Y direction.
                        T_DataTemp = [T_DataTemp;Row,Row,-(BorderMinZ+BorderMaxZ)]; % Tissue domain boundaries in the Z direction.
                        % Establish all interactions with boundaries.
                    else
                        D(Row) = D(Row) + (BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Tissue domain boundaries in the X direction.
                        D(Row) = D(Row) + (BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Tissue domain boundaries in the Y direction.
                        D(Row) = D(Row) + (BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Tissue domain boundaries in the Z direction.
                        T_DataTemp = [T_DataTemp;Row,Row,(BorderMinX+BorderMaxX)*Ax*-H_Out]; % Tissue domain boundaries in the X direction.
                        T_DataTemp = [T_DataTemp;Row,Row,(BorderMinY+BorderMaxY)*Ay*-H_Out]; % Tissue domain boundaries in the Y direction.
                        T_DataTemp = [T_DataTemp;Row,Row,(BorderMinZ+BorderMaxZ)*Az*-H_Out]; % Tissue domain boundaries in the Z direction.
                    % Establish all interactions with boundaries.
                    end % If boundary heat transfer is finite do heat transfer interactions between voxels and boundaries.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Conductive Interactions %%%%%%%%
                    if ~(I==1 || ~DomTot(I-1,J,K)) % If no boundary between current voxel and I-1.
                        T_DataTemp = [T_DataTemp;Row,RowDI,(Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    end
                    % Conduction with I-1.
                    
                    if ~(I==size(DomTot,1) || ~DomTot(I+1,J,K)) % If no boundary between current voxel and I+1.
                        T_DataTemp = [T_DataTemp;Row,RowUI,(Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    end
                    % Conduction with I+1.
                    
                    if ~(J==1 || ~DomTot(I,J-1,K)) % If no boundary between current voxel and J-1.
                        T_DataTemp = [T_DataTemp;Row,RowDJ,(Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    end
                    % Conduction with J-1.
                    
                    if ~(J==size(DomTot,2) || ~DomTot(I,J+1,K)) % If no boundary between current voxel and J+1.
                        T_DataTemp = [T_DataTemp;Row,RowUJ,(Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    end
                    % Conduction with J+1.
                    
                    if ~(K==1 || ~DomTot(I,J,K-1)) % If no boundary between current voxel and K-1.
                        T_DataTemp = [T_DataTemp;Row,RowDK,(Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    end
                    % Conduction with K-1.
                    
                    if ~(K==size(DomTot,3) || ~DomTot(I,J,K+1)) % If no boundary between current voxel and K+1.
                        T_DataTemp = [T_DataTemp;Row,RowUK,(Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    end
                    % Conduction with K+1.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Non-Boundary Interactions %%%%%%%%%%%%%%%%
                else % If voxel does not lie on domain boundary.
                    
                    
                    %%%%%% Establish Conductive Interactions %%%%%%%%
                    T_DataTemp = [T_DataTemp;Row,RowUI,(Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (I+1).
                    T_DataTemp = [T_DataTemp;Row,RowDI,(Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (I-1).
                    T_DataTemp = [T_DataTemp;Row,Row,-(2*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    % Conduction with I-1 and I+1.
                    
                    T_DataTemp = [T_DataTemp;Row,RowUJ,(Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (J+1).
                    T_DataTemp = [T_DataTemp;Row,RowDJ,(Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (J-1).
                    T_DataTemp = [T_DataTemp;Row,Row,-(2*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    % Conduction with J-1 and J+1.
                    
                    T_DataTemp = [T_DataTemp;Row,RowUK,(Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (K+1).
                    T_DataTemp = [T_DataTemp;Row,RowDK,(Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (K-1).
                    T_DataTemp = [T_DataTemp;Row,Row,-(2*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    % Conduction with K-1 and K+1.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Metabolism %%%%%%%%%%%%%%%%%%%%%
                D(Row) = D(Row) - Vol*Q_t; % Tissue domain source for source.
                % Establishes the metabolic heat source term.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Pennes Perfusion %%%%%%%%%%%%%%%
                D(Row) = D(Row) - Vol*Perfusion(I,J,K)*Cp_b*BloodTemp; % Tissue domain interaction for source.
                T_DataTemp = [T_DataTemp;Row,Row,-Vol*Perfusion(I,J,K)*Cp_b]; % Tissue domain interaction for current voxel.
                % Establishes the Pennes Perfusion source term for tissue.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
                T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
                T_DataTemp = []; % Clear T_DataTemp
                % Compiling interactions into T_Super list.
                
            end
        end
    end
end

toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Sparse Matrix %%%%%%%%%%%%%%%%%%%
T_Boolean = logical(T_Super(:,1));
T_Super = T_Super(T_Boolean,:); % Delete any excess rows from T_Super.
T_Solve = sparse(T_Super(:,1),T_Super(:,2),T_Super(:,3)); % Convert P_Super into sparse matrix.
clearvars T_Super T_Boolean % Delete variables to free up memory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
