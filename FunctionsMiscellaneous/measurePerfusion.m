MassFlow = zeros(size(GM_WM));
MassFlowCheck = zeros(size(GM_WM));
for I=1:size(GM_WM,1)
    for J=1:size(GM_WM,2)
        for K=1:size(GM_WM,3)
            if GM_WM(I,J,K)
                    MassFlow(I,J,K) = MassFlow(I,J,K) + abs(Mdot1(I,J,K));

                    if U(I,J,K)>0, MassFlow(I,J,K) = MassFlow(I,J,K) + abs(U(I,J,K))*VoxelSize^2*Rho_b; end
                    if U(I+1,J,K)<0, MassFlow(I,J,K) = MassFlow(I,J,K) + abs(U(I+1,J,K))*VoxelSize^2*Rho_b; end
                    if V(I,J,K)>0, MassFlow(I,J,K) = MassFlow(I,J,K) + abs(V(I,J,K))*VoxelSize^2*Rho_b; end
                    if V(I,J+1,K)<0, MassFlow(I,J,K) = MassFlow(I,J,K) + abs(V(I,J+1,K))*VoxelSize^2*Rho_b; end
                    if W(I,J,K)>0, MassFlow(I,J,K) = MassFlow(I,J,K) + abs(W(I,J,K))*VoxelSize^2*Rho_b; end
                    if W(I,J,K+1)<0, MassFlow(I,J,K) = MassFlow(I,J,K) + abs(W(I,J,K+1))*VoxelSize^2*Rho_b; end
                    
                    MassFlowCheck(I,J,K) = VoxelSize^2*Rho_b*(U(I,J,K)-U(I+1,J,K)+V(I,J,K)-V(I,J+1,K)+W(I,J,K)-W(I,J,K+1)) + Mdot1(I,J,K) + Mdot2(I,J,K);
                
            end
        end
    end
end

MassFlow(~GM_WM) = NaN; % kg/s
MeasuredPerfusion = MassFlow/VoxelSize^3./Rho/Rho_b*10^6/10*60; % ml/100g/min

