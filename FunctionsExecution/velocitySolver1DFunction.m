function [MassFlow, Split] = velocitySolver1DFunction(Vessel,VoxelSize,S,Aavg,L,Rho_b,Visc_b,InletPoints,InletFlows,OutletPoints,OutletFlows)

% function outputs:
%               mass flow rate at each segment (with flow after removing or
%               before adding leakage depending)
% function inputs:
%               Vessel tree structure
%               source list (saved at point but corresponds to branch connection)
%               inlet point list (can be empty, 2 columns [points massflow])
%               outlet point list(can be empty), 2 columns [points pressure], pressure 0 by default)

%%%%%%% Finding Branch Terminations %%%%%%%%%%%%%
BranchTermination = ones(size(Vessel,1),1);
for N=2:numel(BranchTermination)
    BranchTermination(Vessel(N,7)) = 0;
    % The node that the current node connects to cannot be a branch
    % termination.
    if (~isempty(InletPoints) && ismember(N,InletPoints))
        BranchTermination(N) = 0;
        % If point is an inlet then it is not a branch termination.
    end
    if (~isempty(OutletPoints) && ismember(N,OutletPoints))
        BranchTermination(N) = 0;
        % If point is an outlet then it is not a branch termination.
    end
end
% This runs through the tree to find all terminations. For every point it
% checks the defined connection. That connection cannot be a branch
% termination as it has a point connected to it. Therefore all branch
% terminations remain after a single run through of the tree.

if (isempty(InletPoints)||~ismember(1,InletPoints)) && (isempty(OutletPoints)||~ismember(1,OutletPoints)) && numel(find(Vessel(:,7) == 1),1)==1
    BranchTermination(1) = 1;
end
% Check if point 1 lies on a branch termination, inlet or outlet as it
% doesn't have a defined connection. If it is not an inlet, not an outlet
% and has only one connecting node then it is a branch termination.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Set Up Matrix for Pressure Solving %%%%%%
NoPts = size(Vessel,1);
MatrixRows = zeros(10*2*NoPts,3);
RowCount = 1;
D = zeros(2*NoPts,1);

for N=1:NoPts;
    
    if N~=1 % If not node 1. Node 1 has to be dealt with sepeartely
        
        Conn = Vessel(N,7);
        % Find connecting node.
        
        G1 = Rho_b*Aavg(N)^2/(8*pi*Visc_b*VoxelSize*L(N)/2);
        G2 = Rho_b*Aavg(N)^2/(8*pi*Visc_b*VoxelSize*L(N)/2);
        % Create conductance for two halves of the segment
        
        if BranchTermination(N) ~= 1 % If current node is not a branch termination
            
            MatrixRows(RowCount,:) = [2*N,2*N,G2]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N,2*N-1,-G2]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N-1,2*N-1,G2]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N-1,2*N,-G2]; RowCount=RowCount+1;
            % Does interaction with current node and central point
            
            if Conn~=1 || BranchTermination(1) == 0 % If connects to point 1 and point 1 is not a branch termination.
                MatrixRows(RowCount,:) = [2*Conn,2*Conn,G1]; RowCount=RowCount+1;
                MatrixRows(RowCount,:) = [2*Conn,2*N-1,-G1]; RowCount=RowCount+1;
                MatrixRows(RowCount,:) = [2*N-1,2*N-1,G1]; RowCount=RowCount+1;
                MatrixRows(RowCount,:) = [2*N-1,2*Conn,-G1]; RowCount=RowCount+1;
                % Does interaction with connecting node and central point
            else % If connects to point 1 and point 1 is a branch termination.
                MatrixRows(RowCount,:) = [2*Conn,2*Conn,1]; RowCount=RowCount+1;
                MatrixRows(RowCount,:) = [2*Conn,2*N-1,-1]; RowCount=RowCount+1;
                MatrixRows(RowCount,:) = [2*N-1,2*N-1,1]; RowCount=RowCount+1;
                MatrixRows(RowCount,:) = [2*N-1,2*Conn,-1]; RowCount=RowCount+1;
                % Sets central node equal to point 1.There should be no 
                % flow at branch termination.
            end
            
            D(2*N-1) = D(2*N-1) -S(N); % Source term for central point
            D(2*N) = 0; % Nodes have no domain transfer.
            
        else % If current node is a branch termination
            
            MatrixRows(RowCount,:) = [2*N,2*N,1]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N,2*N-1,-1]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N-1,2*N-1,1]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N-1,2*N,-1]; RowCount=RowCount+1;
            % Sets current node equal to central point. There should be no
            % flow at a branch termination.
            
            MatrixRows(RowCount,:) = [2*Conn,2*Conn,G1]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*Conn,2*N-1,-G1]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N-1,2*N-1,G1]; RowCount=RowCount+1;
            MatrixRows(RowCount,:) = [2*N-1,2*Conn,-G1]; RowCount=RowCount+1; 
            % Does interaction with connecting node and central point
            
            D(2*N-1) = D(2*N-1) -S(N); % Source term for central point
            D(2*N) = 0; % Nodes have no domain transfer.
        end
        
    else % If this is node 1.
        MatrixRows(RowCount,:) = [2*N-1,2*N-1,1]; RowCount=RowCount+1;
        MatrixRows(RowCount,:) = [2*N-1,2*N,-1]; RowCount=RowCount+1;
        MatrixRows(RowCount,:) = [2*N,2*N,1]; RowCount=RowCount+1;
        MatrixRows(RowCount,:) = [2*N,2*N-1,-1]; RowCount=RowCount+1;
        % Set node one equal to the theoretical central point on the
        % segment.
        
        D(2*N) = 0; % Nodes have no domain transfer.
    end
    
    if ~isempty(InletPoints)
        if ismember(N,InletPoints)
            D(2*N) = D(2*N)+ InletFlows(InletPoints==N);
        end
    end
    % If there are inlet nodes, set the flowrate at those nodes.
    
    if ~isempty(OutletPoints)
        if ismember(N,OutletPoints)
            D(2*N) = D(2*N) + OutletFlows(OutletPoints==N);
        end
    end
    % If there are outlet nodes, set the flowrate at those nodes.
    
end

N = find(BranchTermination==1,1,'first');
MatrixRows(RowCount,:) = [2*N,2*N,1]; RowCount=RowCount+1;
D(2*N) = 0;
% For a single branch termination, set the pressure to zero to bound the
% solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Solving Linear System %%%%%%%%%%%%%%%%%%%
Matrix_rows_logical = logical(MatrixRows(:,1));
MatrixRows = MatrixRows(Matrix_rows_logical,:);
Matrix = sparse(MatrixRows(:,1),MatrixRows(:,2),MatrixRows(:,3));
P = Matrix\D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%% Converting Pressures into Mass Flows %%%%
MassFlow = zeros(NoPts,2);
Split = zeros(NoPts,2);
% Find Mass Flowrate

MassFlow(1,:) = [0 G2*(P(2*1-1)-P(1*1))];
for N = 2:NoPts
    Conn = Vessel(N,7);

    G1 = Rho_b*Aavg(N)^2/(8*pi*Visc_b*VoxelSize*L(N)/2);
    G2 = Rho_b*Aavg(N)^2/(8*pi*Visc_b*VoxelSize*L(N)/2);
    % Create conductance for two halves of the segment

    MassFlow(N,:) = [G1*(P(2*Conn)-P(2*N-1)) G2*(P(2*N-1)-P(2*N))];
    % Determines mass flowrates for the segment. The first column refers to
    % the flowrate closer to the connecting node. The second column refers
    % to the flowrate closer to the current node (node N). Positive flows
    % directed towards the current node and negative flows are directed
    % towards the connecting node.
    
    if (P(2*Conn) < P(2*N-1) && P(2*N) < P(2*N-1)) || (P(2*Conn) > P(2*N-1) && P(2*N) > P(2*N-1))
        Split(N,:) = [1, G2*abs(P(2*N-1)-P(2*N))/abs(S(N))];
        % Split in flow occurs. Domain mass transfer comes from or flows 
        % towards both sides of the segment. Needs to be checked for proper
        % transfers in temperature solver. The fraction in the second
        % column refers to the fraction of mass transfer that comes from or
        % flows to the current node (node N).
    end
end

