

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


%%%%%% Check if 2 Phases in Porous Domain %%%%%%%
RowAdjustment = 0; % Value that adds rows into the matrix if required. (RowAdjustment = 0 adds no rows).
if Beta23 ~= Inf % If the inter-domain heat transfer is finite.
    RowAdjustment = NumDomTot; % Add rows to solve the tissue and blood phases separately.
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Establishing Number of Rows in Domains %%%
NumDomRows = NumDomTot + RowAdjustment; % Total rows assigned to voxels.
VesselRow = size(Vessel1,1); % Total rows assigned to arteries.
% In the sparse matrix, the rows are divided up by domiain in the following
% order: [tissue -> blood -> arteries -> veins]. These values allow rapid
% identification of rows for arteries and veins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initialising Storage Matrices %%%%%%%%%%%%
T_DataTemp = [];
T_DataRowcount = 1;
T_Super = zeros(1e8,3);
D = zeros(NumDomRows+size(Vessel1,1)+size(Vessel2,1),1);
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
                Row = DomTotConvert(I,J,K); BRow = Row + RowAdjustment;
                if I<size(GM_WM,1), RowUI = DomTotConvert(I+1,J,K); BRowUI = RowUI+RowAdjustment; end
                if I>1,             RowDI = DomTotConvert(I-1,J,K); BRowDI = RowDI+RowAdjustment; end
                if J<size(GM_WM,2), RowUJ = DomTotConvert(I,J+1,K); BRowUJ = RowUJ+RowAdjustment; end
                if J>1,             RowDJ = DomTotConvert(I,J-1,K); BRowDJ = RowDJ+RowAdjustment; end
                if K<size(GM_WM,3), RowUK = DomTotConvert(I,J,K+1); BRowUK = RowUK+RowAdjustment; end
                if K>1,             RowDK = DomTotConvert(I,J,K-1); BRowDK = RowDK+RowAdjustment; end
                % Establishes the row values in the matrix for the current
                % voxel and all neighbouring voxels. If there is only a
                % single phase being solved (because Beta23 = Inf) then all
                % BRow values equal Row values.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Preallocating Searchable Parameters %%%%%%
                Vol=VoxelSize^3;
                Ax=VoxelSize^2;
                Ay=VoxelSize^2;
                Az=VoxelSize^2;
                % Volume and areas of voxels.
                
                if Option_CounterCurrentFlow
                    Porosity_b = 2*Porosity(I,J,K);
                    Porosity_t = 1-2*Porosity_b;
                else
                    Porosity_b = Porosity(I,J,K);
                    Porosity_t = 1-Porosity_b;
                end
                
                if strcmp(TrialOverride,'tracer')
                    Rho_t = 1;
                else
                    Rho_t = Rho(I,J,K);
                end
                Cp_t = Cp(I,J,K);
                Kc_t = Kc(I,J,K);
                Q_t = Q(I,J,K);
                % Establishing searchable parameters. This prevents
                % constantly having to look up values in matrix.
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
                    
                    if Option_AdiabaticBase && K == 1, BorderMinZ=0; end
                    % Set base of model to adiabatic.
                    
                    if Option_AdiabaticArmEnds && K == find(any(any(GM_WM))==1,1,'first'), BorderMinZ=0; end
                    if Option_AdiabaticArmEnds && K == find(any(any(GM_WM))==1,1,'last'), BorderMaxZ=0; end
                    % Set ends of arm model to adiabatic.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Boundary Interactions %%%%%%%%%%
                    if H_Out == Inf % If boundary heat transfer is infinite, set voxel temperature equal to boundary temperature.
                        D(Row) = D(Row) + -(BorderMinX+BorderMaxX)*T_Out; % Tissue domain boundaries in the X direction.
                        D(Row) = D(Row) + -(BorderMinY+BorderMaxY)*T_Out; % Tissue domain boundaries in the Y direction.
                        D(Row) = D(Row) + -(BorderMinZ+BorderMaxZ)*T_Out; % Tissue domain boundaries in the Z direction.
                        D(BRow) = D(BRow) + -(BorderMinX+BorderMaxX)*T_Out; % Blood domain boundaries in the X direction.
                        D(BRow) = D(BRow) + -(BorderMinY+BorderMaxY)*T_Out; % Blood domain boundaries in the Y  direction.
                        D(BRow) = D(BRow) + -(BorderMinZ+BorderMaxZ)*T_Out; % Blood domain boundaries in the Z direction.
                        T_DataTemp = [T_DataTemp;Row,Row,-(BorderMinX+BorderMaxX)]; % Tissue domain boundaries in the X direction.
                        T_DataTemp = [T_DataTemp;Row,Row,-(BorderMinY+BorderMaxY)]; % Tissue domain boundaries in the Y direction.
                        T_DataTemp = [T_DataTemp;Row,Row,-(BorderMinZ+BorderMaxZ)]; % Tissue domain boundaries in the Z direction.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(BorderMinX+BorderMaxX)]; % Blood domain boundaries in the X direction.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(BorderMinY+BorderMaxY)]; % Blood domain boundaries in the Y direction.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(BorderMinZ+BorderMaxZ)]; % Blood domain boundaries in the Z direction.
                        % Establish all interactions with boundaries if
                        % required.
                        
                    else % If boundary heat transfer is finite do heat transfer interactions between voxels and boundaries.
                        D(Row) = D(Row) + Porosity_t*(BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Tissue domain boundaries in the X direction.
                        D(Row) = D(Row) + Porosity_t*(BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Tissue domain boundaries in the Y direction.
                        D(Row) = D(Row) + Porosity_t*(BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Tissue domain boundaries in the Z direction.
                        D(BRow) = D(BRow) + Porosity_b*(BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Blood domain boundaries in the X direction.
                        D(BRow) = D(BRow) + Porosity_b*(BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Blood domain boundaries in the Y  direction.
                        D(BRow) = D(BRow) + Porosity_b*(BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Blood domain boundaries in the Z direction.                      
                        T_DataTemp = [T_DataTemp;Row,Row,Porosity_t*(BorderMinX+BorderMaxX)*Ax*-H_Out]; % Tissue domain boundaries in the X direction.
                        T_DataTemp = [T_DataTemp;Row,Row,Porosity_t*(BorderMinY+BorderMaxY)*Ay*-H_Out]; % Tissue domain boundaries in the Y direction.
                        T_DataTemp = [T_DataTemp;Row,Row,Porosity_t*(BorderMinZ+BorderMaxZ)*Az*-H_Out]; % Tissue domain boundaries in the Z direction.
                        T_DataTemp = [T_DataTemp;BRow,BRow,Porosity_b*(BorderMinX+BorderMaxX)*Ax*-H_Out]; % Blood domain boundaries in the X direction.
                        T_DataTemp = [T_DataTemp;BRow,BRow,Porosity_b*(BorderMinY+BorderMaxY)*Ay*-H_Out]; % Blood domain boundaries in the Y direction.
                        T_DataTemp = [T_DataTemp;BRow,BRow,Porosity_b*(BorderMinZ+BorderMaxZ)*Az*-H_Out]; % Blood domain boundaries in the Z direction.
                        % Establish all interactions with boundaries if
                        % required.
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Conductive Interactions %%%%%%%%
                    if ~(I==1 || ~DomTot(I-1,J,K)) % If no boundary between current voxel and I-1.
                        T_DataTemp = [T_DataTemp;Row,RowDI,(Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTemp = [T_DataTemp;BRow,BRowDI,(Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with I-1 if not a boundary.
                    
                    if ~(I==size(DomTot,1) || ~DomTot(I+1,J,K)) % If no boundary between current voxel and I+1.
                        T_DataTemp = [T_DataTemp;Row,RowUI,(Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTemp = [T_DataTemp;BRow,BRowUI,(Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with I+1 if not a boundary.
                    
                    if ~(J==1 || ~DomTot(I,J-1,K)) % If no boundary between current voxel and J-1.
                        T_DataTemp = [T_DataTemp;Row,RowDJ,(Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTemp = [T_DataTemp;BRow,BRowDJ,(Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with J-1 if not a boundary.
                    
                    if ~(J==size(DomTot,2) || ~DomTot(I,J+1,K)) % If no boundary between current voxel and J+1.
                        T_DataTemp = [T_DataTemp;Row,RowUJ,(Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTemp = [T_DataTemp;BRow,BRowUJ,(Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with J+1 if not a boundary.
                    
                    if ~(K==1 || ~DomTot(I,J,K-1)) % If no boundary between current voxel and K-1.
                        T_DataTemp = [T_DataTemp;Row,RowDK,(Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTemp = [T_DataTemp;BRow,BRowDK,(Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with K-1 if not a boundary.
                    
                    if ~(K==size(DomTot,3) || ~DomTot(I,J,K+1)) % If no boundary between current voxel and K+1.
                        T_DataTemp = [T_DataTemp;Row,RowUK,(Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;Row,Row,-(Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTemp = [T_DataTemp;BRow,BRowUK,(Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTemp = [T_DataTemp;BRow,BRow,-(Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with K+1 if not a boundary.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                %%%%%% Non-Boundary Interactions %%%%%%%%%%%%%%%%
                else % If voxel does not lie on domain boundary.
                    
                    
                    %%%%%% Establish Conductive Interactions %%%%%%%%
                    T_DataTemp = [T_DataTemp;Row,RowUI,(Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (I+1).
                    T_DataTemp = [T_DataTemp;Row,RowDI,(Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (I-1).
                    T_DataTemp = [T_DataTemp;Row,Row,-(2*Porosity_t*Ax*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    T_DataTemp = [T_DataTemp;BRow,BRowUI,(Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel (I+1).
                    T_DataTemp = [T_DataTemp;BRow,BRowDI,(Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel (I-1).
                    T_DataTemp = [T_DataTemp;BRow,BRow,-(2*Porosity_b*Ax*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    % Conduction with I-1 and I+1.
                    
                    T_DataTemp = [T_DataTemp;Row,RowUJ,(Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (J+1).
                    T_DataTemp = [T_DataTemp;Row,RowDJ,(Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (J-1).
                    T_DataTemp = [T_DataTemp;Row,Row,-(2*Porosity_t*Ay*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    T_DataTemp = [T_DataTemp;BRow,BRowUJ,(Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel (J+1).
                    T_DataTemp = [T_DataTemp;BRow,BRowDJ,(Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel (J-1).
                    T_DataTemp = [T_DataTemp;BRow,BRow,-(2*Porosity_b*Ay*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    % Conduction with J-1 and J+1.
                    
                    T_DataTemp = [T_DataTemp;Row,RowUK,(Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (K+1).
                    T_DataTemp = [T_DataTemp;Row,RowDK,(Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for neighbouring voxel (K-1).
                    T_DataTemp = [T_DataTemp;Row,Row,-(2*Porosity_t*Az*(Kc_t/VoxelSize))]; % Tissue domain interaction for current voxel.
                    T_DataTemp = [T_DataTemp;BRow,BRowUK,(Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel (K+1).
                    T_DataTemp = [T_DataTemp;BRow,BRowDK,(Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for neighbouring voxel (K-1).
                    T_DataTemp = [T_DataTemp;BRow,BRow,-(2*Porosity_b*Az*(Kc_b/VoxelSize))]; % Blood domain interaction for current voxel.
                    % Conduction with K-1 and K+1.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% 3D Advection Interactions %%%%%%%%%%%%%%%%
                if U(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRowDI,Ax*Rho_b*Cp_b*abs(U(I,J,K))]; end % Blood domain interaction for upstream voxel.
                if U(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ax*Rho_b*Cp_b*abs(U(I,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from I-1.
                
                if U(I+1,J,K)<0, T_DataTemp = [T_DataTemp;BRow,BRowUI,Ax*Rho_b*Cp_b*abs(U(I+1,J,K))]; end % Blood domain interaction for upstream voxel .
                if U(I+1,J,K)<0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ax*Rho_b*Cp_b*abs(U(I+1,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from I+1.
                
                if V(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRowDJ,Ay*Rho_b*Cp_b*abs(V(I,J,K))]; end % Blood domain interaction for upstream voxel.
                if V(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ay*Rho_b*Cp_b*abs(V(I,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from J-1.
                
                if V(I,J+1,K)<0, T_DataTemp = [T_DataTemp;BRow,BRowUJ,Ay*Rho_b*Cp_b*abs(V(I,J+1,K))]; end % Blood domain interaction for upstream voxel.
                if V(I,J+1,K)<0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ay*Rho_b*Cp_b*abs(V(I,J+1,K))]; end  % Blood domain interaction for current voxel.
                % Advection from J+1.
                
                if W(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRowDK,Az*Rho_b*Cp_b*abs(W(I,J,K))]; end % Blood domain interaction for upstream voxel.
                if W(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRow,-Az*Rho_b*Cp_b*abs(W(I,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from K-1.
                
                if W(I,J,K+1)<0, T_DataTemp = [T_DataTemp;BRow,BRowUK,Az*Rho_b*Cp_b*abs(W(I,J,K+1))]; end % Blood domain interaction for upstream voxel.
                if W(I,J,K+1)<0, T_DataTemp = [T_DataTemp;BRow,BRow,-Az*Rho_b*Cp_b*abs(W(I,J,K+1))]; end  % Blood domain interaction for current voxel.
                % Advection from K+1.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% 3D Advection Interactions (C-C) %%%%%%%%%%
                if Option_CounterCurrentFlow || Option_PseudoCounterCurrentFlow
                    if U2(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRowDI,Ax*Rho_b*Cp_b*abs(U2(I,J,K))]; end % Blood domain interaction for upstream voxel.
                    if U2(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ax*Rho_b*Cp_b*abs(U2(I,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from I-1.
                
                    if U2(I+1,J,K)<0, T_DataTemp = [T_DataTemp;BRow,BRowUI,Ax*Rho_b*Cp_b*abs(U2(I+1,J,K))]; end % Blood domain interaction for upstream voxel .
                    if U2(I+1,J,K)<0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ax*Rho_b*Cp_b*abs(U2(I+1,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from I+1.
                
                    if V2(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRowDJ,Ay*Rho_b*Cp_b*abs(V2(I,J,K))]; end % Blood domain interaction for upstream voxel.
                    if V2(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ay*Rho_b*Cp_b*abs(V2(I,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from J-1.
                
                    if V2(I,J+1,K)<0, T_DataTemp = [T_DataTemp;BRow,BRowUJ,Ay*Rho_b*Cp_b*abs(V2(I,J+1,K))]; end % Blood domain interaction for upstream voxel.
                    if V2(I,J+1,K)<0, T_DataTemp = [T_DataTemp;BRow,BRow,-Ay*Rho_b*Cp_b*abs(V2(I,J+1,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from J+1.
                
                    if W2(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRowDK,Az*Rho_b*Cp_b*abs(W2(I,J,K))]; end % Blood domain interaction for upstream voxel.
                    if W2(I,J,K)>0, T_DataTemp = [T_DataTemp;BRow,BRow,-Az*Rho_b*Cp_b*abs(W2(I,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from K-1.
                
                    if W2(I,J,K+1)<0, T_DataTemp = [T_DataTemp;BRow,BRowUK,Az*Rho_b*Cp_b*abs(W2(I,J,K+1))]; end % Blood domain interaction for upstream voxel.
                    if W2(I,J,K+1)<0, T_DataTemp = [T_DataTemp;BRow,BRow,-Az*Rho_b*Cp_b*abs(W2(I,J,K+1))]; end  % Blood domain interaction for current voxel.
                    % Advection from K+1.
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Metabolism %%%%%%%%%%%%%%%%%%%%%
                D(Row) = D(Row) - Vol*Q_t; % Tissue domain source for source.
                % Establishes the metabolic heat source term.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Pennes Perfusion %%%%%%%%%%%%%%%
                if ~GM_WM(I,J,K) % If outside of the brain domain.
                    D(Row) = D(Row) - Vol*Perfusion(I,J,K)*Cp_b*BloodTemp; % Tissue domain interaction for source.
                    T_DataTemp = [T_DataTemp;Row,Row,-Vol*Perfusion(I,J,K)*Cp_b]; % Tissue domain interaction for current voxel.
                end
                % Establishes the Pennes Perfusion source term for tissue
                % outside of the brain domain.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% 3D Inter-Domain Heat Transfer %%%%%%%%%%%%
                if ~isinf(Beta23) % If the heat transfer coefficient is finite.
                    if Porosity_b ~= 0 % If there are two phases present.
                        T_DataTemp = [T_DataTemp;Row,Row,-Beta23]; % Tissue domain interaction for current voxel (tissue domain).
                        T_DataTemp = [T_DataTemp;Row,BRow,Beta23]; % Tissue domain interaction for current voxel (blood domain).
                        T_DataTemp = [T_DataTemp;BRow,Row,Beta23]; % Blood domain interaction for current voxel (tissue domain).
                        T_DataTemp = [T_DataTemp;BRow,BRow,-Beta23]; % Blood domain interaction for current voxel (blood domain).
                        % Establishes the inter-domain heat transfer
                        % between the tissue phase and the blood phase.
                        
                    else % If only tissue phase is present.
                        T_DataTemp = [T_DataTemp;BRow,Row,2]; % Blood domain interaction for current voxel (tissue domain).
                        T_DataTemp = [T_DataTemp;BRow,BRow,-2]; % Blood domain interaction for current voxel (blood domain).
                        % This sets the blood phase equal to the tissue
                        % phase. This avoids any voxels remaining undefined
                        % if the porosity is 0 (eg outside the brain
                        % domain).
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Transient Heat Transfer %%%%%%%%%%%%%%%%%%
                if Option_TransientSolve
                    T_DataTemp = [T_DataTemp;Row,Row,-2*Vol*Porosity_t*Rho_t*Cp_t/Timestep];
                    T_DataTemp = [T_DataTemp;BRow,BRow,-2*Vol*Porosity_b*Rho_b*Cp_b/Timestep];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Artery Voxel Interactions %%%%%%%%%%%%%%%%
                if ~isempty(Volume2Vessel1{I,J,K}) % If this voxel has intersecting artery segments.
                    for mn = 1:size(Volume2Vessel1{I,J,K},1) % For every intersecting artery segment.
                        
                        %%%%%% Preallocating Searchable Parameters %%%%%%
                        Node = Volume2Vessel1{I,J,K}(mn,1); % Define arterial node.
                        NodeMdot = Volume2Vessel1{I,J,K}(mn,2); % Define inter-domain mass transfer from node segment
                        Beta1_t = Beta1(Node); % Define inter-domain heat transfer coefficient from node segment.
                        % Establishing searchable parameters. This prevents
                        % constantly having to look up values in matrix.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        %%%%%% Split Condition  %%%%%%%%%%%%%%%%%%%%%%%%%
                        if SplitArt(Node,1) == 1
                            % If flow is split along this segment, blood
                            % arrives from both nodes of the segment and no
                            % heat transfer occurs.
                            
                            NodeConn = Vessel1(Node,7); % Define connecting node.
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+Node,SplitArt(Node,2)*NodeMdot*Cp_b]; % Artery domain interaction for node.
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+NodeConn,(1-SplitArt(Node,2))*NodeMdot*Cp_b]; % Artery domain interaction for connecting node.
                            T_DataTemp = [T_DataTemp;BRow,BRow,-NodeMdot*Cp_b]; % Blood domain interaction for current voxel.
                            % Inter-domain convection from mass transfer if
                            % there is a split on the segment. Here blood
                            % is taken porportionally from both nodes in
                            % the segment. No further heat transfer occurs
                            % in the segment.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                        %%%%%% Non-Split Condition  %%%%%%%%%%%%%%%%%%%%%
                        else % If no split on segment.
                            
                            
                            %%%%%% Determine Direction of Flow  %%%%%%%%%%%%%
                            if FdotArt(Node,2) < 0
                                NodeDown = Vessel1(Node,7); % Downstream node.
                                NodeUp = Node; % Upstream node.
                            else
                                NodeDown = Node; % Downstream node.
                                NodeUp = Vessel1(Node,7); % Upstream node.
                            end
                            % This determines which is the upstream and
                            % downstream nodes of the line segments. NodeUp
                            % is defined as the upstream and NodeDown is
                            % downstream.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Convection %%%%%%%%%%%%%%%%%%
                            % Blood source from leakage of arteries
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+NodeUp,NodeMdot*Cp_b]; % Blood domain interaction for upstream arterial node.
                            T_DataTemp = [T_DataTemp;BRow,BRow,-NodeMdot*Cp_b]; % Blood domain interaction for current voxel.
                            % Inter-domain convection from mass transfer.
                            % This comes from the upstream node only.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Tissue) %%%%%%%%%
                            T_DataTemp = [T_DataTemp;Row,NumDomRows+NodeDown,Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for downstream node.
                            T_DataTemp = [T_DataTemp;Row,Row,-Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTemp = [T_DataTemp;Row,NumDomRows+NodeUp,Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for upstream node.
                            T_DataTemp = [T_DataTemp;Row,Row,-Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,Row,Porosity_t*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,NumDomRows+NodeDown,-Porosity_t*Beta1_t*0.5]; % Downstream node interaction for downstream node.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,Row,Porosity_t*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,NumDomRows+NodeUp,-Porosity_t*Beta1_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % from tissue domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Blood) %%%%%%%%%%
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+NodeDown,Porosity_b*Beta1_t*0.5]; % Blood domain interaction for downstream node.
                            T_DataTemp = [T_DataTemp;BRow,BRow,-Porosity_b*Beta1_t*0.5]; % Blood domain interaction for upstream node.
                            
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+NodeUp,Porosity_b*Beta1_t*0.5]; % Blood domain interaction for upstream node.
                            T_DataTemp = [T_DataTemp;BRow,BRow,-Porosity_b*Beta1_t*0.5]; % Blood domain interaction for current voxel.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,BRow,Porosity_b*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,NumDomRows+NodeDown,-Porosity_b*Beta1_t*0.5]; % Downstream node interaction for downstream node.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,BRow,Porosity_b*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+NodeDown,NumDomRows+NodeUp,-Porosity_b*Beta1_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % from blood domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Vein Voxel Interactions %%%%%%%%%%%%%%%%%%
                if ~isempty(Volume2Vessel2{I,J,K}) % If this voxel has intersecting vein segments.
                    for mn = 1:size(Volume2Vessel2{I,J,K},1) % For every intersecting vein segment.
                        
                        %%%%%% Preallocating Searchable Parameters %%%%%%
                        Node = Volume2Vessel2{I,J,K}(mn,1);
                        Beta2_t = Beta2(Node);
                        % Establishing searchable parameters. This prevents
                        % constantly having to look up values in matrix.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        %%%%%% Non-Split Condition  %%%%%%%%%%%%%%%%%%%%%
                        if SplitVein(Node,1) ~= 1
                            % If flow is split along this segment, no heat
                            % transfer occurs.
                            
                            
                            %%%%%% Determine Direction of Flow  %%%%%%%%%%%%%
                            if FdotVein(Node,2) < 0
                                NodeDown = Vessel2(Node,7);
                                NodeUp = Node;
                            else
                                NodeDown = Node;
                                NodeUp = Vessel2(Node,7);
                            end
                            % This determines which is the upstream and
                            % downstream nodes of the line segments. NodeUp
                            % is defined as the upstream and NodeDown is
                            % downstream.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Tissue) %%%%%%%%%
                            T_DataTemp = [T_DataTemp;Row,NumDomRows+VesselRow+NodeDown,Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for downstream node.
                            T_DataTemp = [T_DataTemp;Row,Row,-Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTemp = [T_DataTemp;Row,NumDomRows+VesselRow+NodeUp,Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for upstream node.
                            T_DataTemp = [T_DataTemp;Row,Row,-Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,Row,Porosity_t*Beta2_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeDown,-Porosity_t*Beta2_t*0.5]; % Downstream node interaction for downstream node.
                            %
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,Row,Porosity_t*Beta2_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeUp,-Porosity_t*Beta2_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % from tissue domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Blood) %%%%%%%%%%
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+VesselRow+NodeDown,Porosity_b*Beta2_t*0.5]; % Blood domain interaction for downstream node.
                            T_DataTemp = [T_DataTemp;BRow,BRow,-Porosity_b*Beta2_t*0.5]; % Blood domain interaction for upstream node.
                            
                            T_DataTemp = [T_DataTemp;BRow,NumDomRows+VesselRow+NodeUp,Porosity_b*Beta2_t*0.5]; % Blood domain interaction for upstream node.
                            T_DataTemp = [T_DataTemp;BRow,BRow,-Porosity_b*Beta2_t*0.5]; % Blood domain interaction for current voxel.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,BRow,Porosity_b*Beta2_t*0.5];  % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeDown,-Porosity_b*Beta2_t*0.5]; % Downstream node interaction for downstream node.
                            
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,BRow,Porosity_b*Beta2_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeUp,-Porosity_b*Beta2_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % Inter-domain heat transfer to the segment
                            % from blood domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
                T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
                T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
                T_DataTemp = []; % Clear T_DataTemp
                % Compiling interactions into T_Super list.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
            end
        end
    end
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Artery Segment Interactions %%%%
disp('Compiling Interactions within Arteries')
tic % Start timing.

for Node = 1:size(Vessel2Volume1,1)
    
    %%%%%% Inlet Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(Node,InletPoints)
        Index = find(Node==InletPoints); % Find which inlet is being referred to.
        T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node,2]; % Set node equal to InletTemp.
        D(NumDomRows+Node) = InletTemp(Index); % Set node equal to InletTemp.
        % If node is an inlet, set temperature to corresponding InletTemp.
        
        if ~SplitArt(Node,1) % If current segment is not split.
            NodeConn = Vessel1(Node,7); % Establish connecting node.
            if FdotArt(Node,1) < 0 % If flow is directed towards connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+Node,abs(FdotArt(Node,1))*Cp_b]; % Connecting node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+NodeConn,-abs(FdotArt(Node,1))*Cp_b]; % Connecting node interaction for connecting node.
            end
        end
        % Do advection flow for the downstream node of current segment.
        % Flow cannot be directed towards an inlet so no upstream advection
        % is performed.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    %%%%%% Non-Inlet Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        
        %%%%%% Arterial Advection Interactions %%%%%%%%%%
        if Node~=1 % Node 1 does not contain a segment.
            if ~SplitArt(Node,1) % If current segment is not split.
                NodeConn = Vessel1(Node,7); % Establish connecting node.
                if FdotArt(Node,2) > 0 % If flow is directed towards current node.
                    T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+NodeConn,abs(FdotArt(Node,2))*Cp_b]; % Current node interaction for connecting node.
                    T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node,-abs(FdotArt(Node,2))*Cp_b]; % Current node interaction for current node.
                end
                if FdotArt(Node,1) < 0 % If flow is directed towards connecting node.
                    T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+Node,abs(FdotArt(Node,1))*Cp_b]; % Connecting node interaction for current node.
                    T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+NodeConn,-abs(FdotArt(Node,1))*Cp_b]; % Connecting node interaction for connecting node.
                end
            end
        end
        % Do advection flow for the downstream node of current segment.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%% Arterial Conduction Interactions %%%%%%%%%
        if Node~=1 % Node 1 does not contain a segment.
            NodeConn = Vessel1(Node,7); % Establish connecting node.
            if ~(NodeConn == 1 && ismember(1,InletPoints)) % Imbalance occurs if there is conduction to an inlet point.
                T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+NodeConn, Aavg1(Node)*(Kc_b/L1(Node))]; % Current node interaction for connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node, -Aavg1(Node)*(Kc_b/L1(Node))]; % Current node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+Node, Aavg1(Node)*(Kc_b/L1(Node))]; % Connecting node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+NodeConn, -Aavg1(Node)*(Kc_b/L1(Node))]; % Connecting node interaction for current node.
            
            end
        end
        % Do conduction for the segment. Both current and connecting nodes
        % are updated here so that connections do not have to be searched
        % for in each node.
        % Conduction should not happen to an inlet point as it is defined
        % directly as the inlet temperature.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%% Transient Heat Transfer %%%%%%%%%%%%%%%%%%
    if Option_TransientSolve
        T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node,-2*Vol1(Node)*Rho_b*Cp_b/Timestep];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
    %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
    if size(T_DataTemp,1) ~= 0
        T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
        T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
        T_DataTemp = []; % Clear T_DataTemp.
        % Compiling interactions into T_Super list.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Venous Segment Interactions %%%%
disp('Compiling Interactions within Veins')
tic % Start timing.

% Venous Vessel Heat Flow
for Node = 1:size(Vessel2Volume2,1)
    
    %%%%%% Venous Advection Interactions %%%%%%%%%%%%
    if Node~=1 % Node 1 does not contain a segment.
        if ~SplitVein(Node,1) % If current segment is not split.
            NodeConn = Vessel2(Node,7); % Establish connecting node.
            if FdotVein(Node,1) > 0 % If flow is directed towards current node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+NodeConn,abs(FdotVein(Node,1))*Cp_b]; % Current node interaction for connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node,-abs(FdotVein(Node,1))*Cp_b]; % Current node interaction for current node.
            end
            if FdotVein(Node,2) < 0 % If flow is directed towards connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+Node,abs(FdotVein(Node,2))*Cp_b]; % Connecting node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+NodeConn,-abs(FdotVein(Node,2))*Cp_b]; % Connecting node interaction for current node.
            end
        end
    end
    % Do advection flow for the downstream node of current segment.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%% Venous Conduction Interactions %%%%%%%%%%%
    if Node ~= 1 % Node 1 does not contain a segment.
        NodeConn = Vessel2(Node,7); % Establish connecting node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+NodeConn, Aavg2(Node)*(Kc_b/L2(Node))]; % Current node interaction for connecting node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node, -Aavg2(Node)*(Kc_b/L2(Node))]; % Current node interaction for current node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+Node, Aavg2(Node)*(Kc_b/L2(Node))]; % Connecting node interaction for current node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+NodeConn, -Aavg2(Node)*(Kc_b/L2(Node))]; % Connecting node interaction for current node.
    end
    % Do conduction for the segment. Both current and connecting nodes
    % are updated here so that connections do not have to be searched
    % for in each node.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Inter-Domain Convection %%%%%%%%%%%%%%%%%%
    if Node~=1 % Node 1 does not contain a segment.
        
        if FdotVein(Node,1) > 0
            NodeDown = Node;
        else
            NodeDown = NodeConn;
        end
        % This determines which is the downstream node (NodeDown) of the
        % line segments.
        
        for NumInts = 1:size(Vessel2Volume2{Node},1) % For all intersections of venous segment.
            I = Vessel2Volume2{Node}(NumInts,1); % Get location I for ntersecting voxel.
            J = Vessel2Volume2{Node}(NumInts,2); % Get location J for ntersecting voxel.
            K = Vessel2Volume2{Node}(NumInts,3); % Get location K for ntersecting voxel.
            BRow = DomTotConvert(I,J,K) + RowAdjustment; % Find corresponding row value for intersecting voxel (blood domain).
            VoxMdot = abs(Vessel2Volume2{Node}(NumInts,4)); % Get mass transfer for ntersecting voxel.
            
            if SplitVein(Node,1) % If segment is split.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,BRow,(SplitVein(Node,2))*VoxMdot*Cp_b]; % Current node interaction for blood domain.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node,-(SplitVein(Node,2))*VoxMdot*Cp_b]; % Current node interaction for current node.
                
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,BRow,(1-SplitVein(Node,2))*VoxMdot*Cp_b]; % Connecting node interaction for blood domain.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+NodeConn,-(1-SplitVein(Node,2))*VoxMdot*Cp_b]; % Connecting node interaction for connecting node.
                % Deliver the corresponding amount of blood to both nodes
                % on the segment.
                
            else % If segment is not split.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,BRow,VoxMdot*Cp_b]; % Downstream node interaction for blood domain.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeDown,-VoxMdot*Cp_b]; % Downstream node interaction for downstream node.
                % Deliver the blood to the downstream node only.
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Transient Heat Transfer %%%%%%%%%%%%%%%%%%
    if Option_TransientSolve
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node,-2*Vol2(Node)*Rho_b*Cp_b/Timestep];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
    if size(T_DataTemp,1) ~= 0
        T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
        T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
        T_DataTemp = []; % Clear T_DataTemp.
        % Compiling interactions into T_Super list.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Sparse Matrix %%%%%%%%%%%%%%%%%%%
T_Boolean = logical(T_Super(:,1));
T_Super = T_Super(T_Boolean,:); % Delete any excess rows from T_Super.
T_Super(:,3) = 0.5*T_Super(:,3);
T_Solve = sparse(T_Super(:,1),T_Super(:,2),T_Super(:,3)); % Convert P_Super into sparse matrix.
clearvars T_Super T_Boolean % Delete variables to free up memory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
