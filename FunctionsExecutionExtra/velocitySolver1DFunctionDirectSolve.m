function [MassFlow, Split] = velocitySolver1DFunctionDirectSolve(Vessel,S,InletPoints,InletFlows,OutletPoints,OutletFlows)

% function outputs:
%               Mass flow rate at each segment (with flow after removing or
%               before adding leakage depending)
%               Condition on whether the segment has split flow or not
% function inputs:
%               Vessel tree structure
%               source list (saved at point but corresponds to branch connection)
%               Inlet point list (can be empty)
%               Inlet flowrates (can be empty)
%               Outlet point list(can be empty)
%               Outlet flowrates (can be empty)

%%%%%%% Set Up Matrix for Flow Solving %%%%%%
NoPts = size(Vessel,1);
MatrixRows = zeros(4*NoPts,3);
RowCount = 1;
D = zeros(2*NoPts,1);


for N=1:NoPts;
    
    Conn = Vessel(N,7); % Connection of segment
    
    MatrixRows(RowCount,:) = [2*N,2*N,-1]; RowCount=RowCount+1;
    
    MatrixRows(RowCount,:) = [2*N-1,2*N-1,-1]; RowCount=RowCount+1;
    MatrixRows(RowCount,:) = [2*N-1,2*N,1]; RowCount=RowCount+1;
    
    if N~=1
        MatrixRows(RowCount,:) = [2*Conn,2*N-1,1]; RowCount=RowCount+1;
    end
    
    D(2*N-1) = D(2*N-1) - S(N); % Source term for central point
    D(2*N) = 0; % Nodes have no domain transfer.

    if ~isempty(InletPoints)
        if ismember(N,InletPoints)
            D(2*N) = D(2*N) + InletFlows(InletPoints==N);
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


%%%%%%% Solving Linear System %%%%%%%%%%%%%%%%%%%
Matrix_rows_logical = logical(MatrixRows(:,1));
MatrixRows = MatrixRows(Matrix_rows_logical,:);
Matrix = sparse(MatrixRows(:,1),MatrixRows(:,2),MatrixRows(:,3));
MatrixSolve = Matrix\D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Extracting Results %%%%%%%%%%%%%%%%%%%%%%
MassFlow = zeros(NoPts,2);
Split = zeros(NoPts,2);
% Find Mass Flowrate

for N = 1:NoPts
    MassFlow(N,:) = [MatrixSolve(2*N-1),MatrixSolve(2*N)];
    if sign(MassFlow(N,1))*sign( MassFlow(N,2)) > 0
        Split(N,:) = [1, abs(MassFlow(N,2))/abs(S(N))];
    end
    
end

end