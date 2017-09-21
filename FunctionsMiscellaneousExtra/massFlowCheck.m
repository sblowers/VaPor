% Arteries

MassFlowArteriesCheck = zeros(2*size(Vessel1,1),1);
for N = 1:size(Vessel1,1)
    
    Conn = Vessel1(N,7);
    
    if N~=1
        if FdotArt(N,1) > 0
            MassFlowArteriesCheck(2*N-1) = MassFlowArteriesCheck(2*N-1) + abs(FdotArt(N,1));
            MassFlowArteriesCheck(2*Conn) = MassFlowArteriesCheck(2*Conn) - abs(FdotArt(N,1));
        else
            MassFlowArteriesCheck(2*N-1) = MassFlowArteriesCheck(2*N-1) - abs(FdotArt(N,1));
            MassFlowArteriesCheck(2*Conn) = MassFlowArteriesCheck(2*Conn) + abs(FdotArt(N,1));
        end
        
        if FdotArt(N,2) > 0
            MassFlowArteriesCheck(2*N) = MassFlowArteriesCheck(2*N) + abs(FdotArt(N,2));
            MassFlowArteriesCheck(2*N-1) = MassFlowArteriesCheck(2*N-1) - abs(FdotArt(N,2));
        else
            MassFlowArteriesCheck(2*N) = MassFlowArteriesCheck(2*N) - abs(FdotArt(N,2));
            MassFlowArteriesCheck(2*N-1) = MassFlowArteriesCheck(2*N-1) + abs(FdotArt(N,2));
        end
        
        MassFlowArteriesCheck(2*N-1) = MassFlowArteriesCheck(2*N-1) - MdotArt(N);
    end
    
    if ismember(N,InletPoints)
        MassFlowArteriesCheck(2*N) = MassFlowArteriesCheck(2*N) + InletFlows(InletPoints==N);
    end
    
end

% Veins

MassFlowVeinsCheck = zeros(2*size(Vessel2,1),1);
for N = 1:size(Vessel2,1)
    
    Conn = Vessel2(N,7);
    
    if N~=1
        if FdotVein(N,1) > 0
            MassFlowVeinsCheck(2*N-1) = MassFlowVeinsCheck(2*N-1) + abs(FdotVein(N,1));
            MassFlowVeinsCheck(2*Conn) = MassFlowVeinsCheck(2*Conn) - abs(FdotVein(N,1));
        else
            MassFlowVeinsCheck(2*N-1) = MassFlowVeinsCheck(2*N-1) - abs(FdotVein(N,1));
            MassFlowVeinsCheck(2*Conn) = MassFlowVeinsCheck(2*Conn) + abs(FdotVein(N,1));
        end
        
        if FdotVein(N,2) > 0
            MassFlowVeinsCheck(2*N) = MassFlowVeinsCheck(2*N) + abs(FdotVein(N,2));
            MassFlowVeinsCheck(2*N-1) = MassFlowVeinsCheck(2*N-1) - abs(FdotVein(N,2));
        else
            MassFlowVeinsCheck(2*N) = MassFlowVeinsCheck(2*N) - abs(FdotVein(N,2));
            MassFlowVeinsCheck(2*N-1) = MassFlowVeinsCheck(2*N-1) + abs(FdotVein(N,2));
        end
        
        MassFlowVeinsCheck(2*N-1) = MassFlowVeinsCheck(2*N-1) - MdotVein(N);
    end
    
    if ismember(N,OutletPoints)
        MassFlowVeinsCheck(2*N) = MassFlowVeinsCheck(2*N) + OutletFlows(OutletPoints==N);
    end
    
end