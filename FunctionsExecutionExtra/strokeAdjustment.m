disp('%%%%%%% Adjusting Flowrates from Stroke %%%%%%%%%%%') % Display
disp('Creating Network Trees')
tic
numArt = size(Vessel1,1);
numVein = size(Vessel2,1);
nxnynz = numel(GM_WM(GM_WM));
row_convert = find(GM_WM); % goes from row number back to GM_WM
row_convert2 = (1:nxnynz) + numArt;
GM_WM_convert = zeros(size(GM_WM));
GM_WM_convert(row_convert) = row_convert2;

error_amt = 1e-15;

Tree = cell(numArt+nxnynz+numVein,1);
% data stored in [location, flowrate to that location]
% for arteries and veins [location, flowrate to that location, vessel segment flow alteration refers to]
% for 3D flow [location, flowrate to that location, number 1 to 6 to reflect the directional flow] 


% find downstream flows and weights for arteries
for art_loop = 1:numArt
    
    if FdotArt(art_loop,1) <= 0 
        % probablility of travelling down artery length
        if abs(FdotArt(art_loop,1)) > 0
            Tree{art_loop}(end+1,:) = [Vessel1(art_loop,7),abs(FdotArt(art_loop,1)),art_loop];
        end
        
        % probablility of transferring domain on artery length
        if ~isempty(Vessel2Volume1{art_loop})
            for n = 1:size(Vessel2Volume1{art_loop},1)
                if Vessel2Volume1{art_loop}(n,4) ~= 0;
                    Tree{art_loop}(end+1,:) = [GM_WM_convert(Vessel2Volume1{art_loop}(n,1),Vessel2Volume1{art_loop}(n,2),Vessel2Volume1{art_loop}(n,3))...
                        ,Vessel2Volume1{art_loop}(n,4),art_loop];
                end
            end
        end
        
    elseif FdotArt(art_loop,2) >= 0
        
        % probablility of travelling down artery length
        if abs(FdotArt(art_loop,2)) > 0 && abs(FdotArt(art_loop,2)) > error_amt
            Tree{Vessel1(art_loop,7)}(end+1,:) = [art_loop,abs(FdotArt(art_loop,2)),art_loop];
        end
        
        % probablility of transferring domain on artery length
        if ~isempty(Vessel2Volume1{art_loop})
            for n = 1:size(Vessel2Volume1{art_loop},1)
                if Vessel2Volume1{art_loop}(n,4) ~= 0;
                    Tree{Vessel1(art_loop,7)}(end+1,:) = [GM_WM_convert(Vessel2Volume1{art_loop}(n,1),Vessel2Volume1{art_loop}(n,2),Vessel2Volume1{art_loop}(n,3))...
                        ,Vessel2Volume1{art_loop}(n,4),art_loop];
                end
            end
        end
        
    end
    
end

% find downstream flows and weights for domain
for domain_loop = 1:nxnynz
    down_temp = [];
    
    [I,J,K] = ind2sub(size(GM_WM),row_convert(domain_loop));
    
    Ax = VoxelSize^2; Ay = VoxelSize^2; Az = VoxelSize^2;
    if U(I,J,K)<0, down_temp(end+1,:) = [GM_WM_convert(I-1,J,K),Ax*Rho_b*abs(U(I,J,K)),1]; end 
    if U(I+1,J,K)>0, down_temp(end+1,:) = [GM_WM_convert(I+1,J,K),Ax*Rho_b*abs(U(I+1,J,K)),2]; end
    if V(I,J,K)<0, down_temp(end+1,:) = [GM_WM_convert(I,J-1,K),Ay*Rho_b*abs(V(I,J,K)),3]; end
    if V(I,J+1,K)>0, down_temp(end+1,:) = [GM_WM_convert(I,J+1,K),Ay*Rho_b*abs(V(I,J+1,K)),4]; end
    if W(I,J,K)<0, down_temp(end+1,:) = [GM_WM_convert(I,J,K-1),Az*Rho_b*abs(W(I,J,K)),5]; end
    if W(I,J,K+1)>0, down_temp(end+1,:) = [GM_WM_convert(I,J,K+1),Az*Rho_b*abs(W(I,J,K+1)),6]; end
   
    if ~isempty(Volume2Vessel2{I,J,K})
        for n = 1:size(Volume2Vessel2{I,J,K},1)
            feed = Volume2Vessel2{I,J,K}(n,1);
            if SplitVein(feed,1)~=0
                down_temp(end+1,:) = [feed+numArt+nxnynz,SplitVein(n,2)*abs(Volume2Vessel2{I,J,K}(n,2)),Volume2Vessel2{I,J,K}(n,1)];
                down_temp(end+1,:) = [Vessel2(feed,7)+numArt+nxnynz,(1-SplitVein(n,2))*abs(Volume2Vessel2{I,J,K}(n,2)),Volume2Vessel2{I,J,K}(n,1)];
            else
                if FdotVein(feed,2) <= 0
                    feed = Vessel2(feed,7);
                end
                down_temp(end+1,:) = [feed+numArt+nxnynz,abs(Volume2Vessel2{I,J,K}(n,2)),Volume2Vessel2{I,J,K}(n,1)];
                
            end
        end
    end
    
    if ~isempty(down_temp)
        Tree{domain_loop+numArt} = [Tree{domain_loop+numArt};down_temp];
    end
    
end

% find downstream flows and weights for veins
for vein_loop = 1:numVein
    
    if FdotVein(vein_loop,2) < 0
        Tree{vein_loop+numArt+nxnynz}(end+1,:) = [Vessel2(vein_loop,7)+numArt+nxnynz,abs(FdotVein(vein_loop,2)),vein_loop];
    elseif FdotVein(vein_loop,1) > 0
        Tree{Vessel2(vein_loop,7)+numArt+nxnynz}(end+1,:) = [vein_loop+numArt+nxnynz,abs(FdotVein(vein_loop,1)),vein_loop];
    end
    
end
toc

% Do the backwards flow for stroke point in arterial domain only
tic
Tree2 = cell(numArt,1);

for art_loop = 1:numArt
    
    if FdotArt(art_loop,1) < 0 && abs(FdotArt(art_loop,1)) > error_amt
%             Tree2{art_loop}(end+1,:) = [Vessel1(art_loop,7),abs(FdotArt(art_loop,1)),art_loop]; 
            Tree2{Vessel1(art_loop,7)}(end+1,:) = [art_loop,abs(FdotArt(art_loop,1)),art_loop]; 
            
    elseif FdotArt(art_loop,2) > 0 && abs(FdotArt(art_loop,2)) > error_amt
%             Tree2{Vessel1(art_loop,7)}(end+1,:) = [art_loop,abs(FdotArt(art_loop,2)),art_loop];
            Tree2{art_loop}(end+1,:) = [Vessel1(art_loop,7),abs(FdotArt(art_loop,2)),art_loop];    
            
    end
    
end
toc


disp('Adjusting Arterial Flowrates (Downstream)')
tic
% start_position = 1;
% start_position = inlet(2,1); % arterial node of stroke occurence
% start_position = 825; % left middle cerebral artery
start_position = StrokeLocation;
Destinations = Tree{start_position};

stroke_source = zeros(size(GM_WM));

while ~isempty(Destinations)
    Destinations_new = [];
    
    for n = 1:size(Destinations,1)
        if Destinations(n,1) <= numArt
            
            % Reducing flow down that stream
            if FdotArt(Destinations(n,3),1) <= 0 && FdotArt(Destinations(n,3),2) <= 0
                FdotArt(Destinations(n,3),:) = FdotArt(Destinations(n,3),:) + Destinations(n,2); % Upstream and Downstream
            else
                FdotArt(Destinations(n,3),:) = FdotArt(Destinations(n,3),:) - Destinations(n,2); % Upstream and Downstream
            end

            if ~isempty(Tree{Destinations(n,1)})
                Destination_alter = [Tree{Destinations(n,1)}(:,1),...
                    (Destinations(n,2)*Tree{Destinations(n,1)}(:,2)/sum(Tree{Destinations(n,1)}(:,2))),...
                    Tree{Destinations(n,1)}(:,3)];
                Destinations_new = [Destinations_new;Destination_alter];
            end
        else
            
            % Reducing flow down that stream
            if FdotArt(Destinations(n,3),1) <= 0 && FdotArt(Destinations(n,3),2) <= 0
                FdotArt(Destinations(n,3),2) = FdotArt(Destinations(n,3),2) + Destinations(n,2); % Upstream only
            else
                FdotArt(Destinations(n,3),1) = FdotArt(Destinations(n,3),1) - Destinations(n,2); % Upstream only
            end
            
            [I,J,K] = ind2sub(size(GM_WM),row_convert(Destinations(n,1)-numArt));
            
            index = find(Volume2Vessel1{I,J,K}(:,1)==Destinations(n,3));
            Volume2Vessel1{I,J,K}(index,2) = Volume2Vessel1{I,J,K}(index,2) - Destinations(n,2);
            index = find(Vessel2Volume1{Destinations(n,3)}(:,1) == I & Vessel2Volume1{Destinations(n,3)}(:,2) == J & Vessel2Volume1{Destinations(n,3)}(:,3) == K);
            Vessel2Volume1{Destinations(n,3)}(index,4) = Vessel2Volume1{Destinations(n,3)}(index,4) - Destinations(n,2);
            
            Mdot1(I,J,K) = Mdot1(I,J,K) - Destinations(n,2);
            
            stroke_source(I,J,K) = stroke_source(I,J,K) + Destinations(n,2);
        end
    end
    Destinations = Destinations_new;
%     size(Destinations,1)
end

toc
disp('Adjusting Porous Flowrates (Downstream)')
tic
stroke_source_veins = zeros(size(Vessel2,1),1);

stroke_source(stroke_source<error_amt) = 0;
stroke_source2 = zeros(size(stroke_source));
while sum(sum(sum(stroke_source))) > 0
    for I = 1:size(stroke_source,1)
        for J = 1:size(stroke_source,2)
            for K = 1:size(stroke_source,3)
                if stroke_source(I,J,K)
                    Destinations = Tree{GM_WM_convert(I,J,K)};
                    outflow = sum(Destinations(:,2));
                    frac = Destinations(:,2)/outflow;
                    for n = 1:size(Destinations,1)
                        if  Destinations(n,1) <= (numArt + nxnynz)

                            [i2,j2,k2] = ind2sub(size(GM_WM),row_convert(Destinations(n,1)-numArt));

                            if Destinations(n,3) == 1, U(I,J,K) = U(I,J,K) + 1/Rho_b*1/VoxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 2, U(I+1,J,K) = U(I+1,J,K) - 1/Rho_b*1/VoxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 3, V(I,J,K) = V(I,J,K) + 1/Rho_b*1/VoxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 4, V(I,J+1,K) = V(I,J+1,K) - 1/Rho_b*1/VoxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 5, W(I,J,K) = W(I,J,K) + 1/Rho_b*1/VoxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 6, W(I,J,K+1) = W(I,J,K+1) - 1/Rho_b*1/VoxelSize^2*frac(n)*stroke_source(I,J,K); end
                            
                            
                            stroke_source2(i2,j2,k2) = stroke_source2(i2,j2,k2) + frac(n)*stroke_source(I,J,K);
                        else
                            
                            if FdotVein(Destinations(n,3),2) > 0
                                FdotVein(Destinations(n,3),2) = FdotVein(Destinations(n,3),2) - frac(n)*stroke_source(I,J,K);
                            else
                                FdotVein(Destinations(n,3),1) = FdotVein(Destinations(n,3),1) + frac(n)*stroke_source(I,J,K);
                            end
                            
                            
                            index = find(Volume2Vessel2{I,J,K}(:,1)==Destinations(n,3));
                            Volume2Vessel2{I,J,K}(index,2) = Volume2Vessel2{I,J,K}(index,2) + frac(n)*stroke_source(I,J,K);
                            index = find(Vessel2Volume2{Destinations(n,3)}(:,1) == I & Vessel2Volume2{Destinations(n,3)}(:,2) == J & Vessel2Volume2{Destinations(n,3)}(:,3) == K);
                            Vessel2Volume2{Destinations(n,3)}(index,4) = Vessel2Volume2{Destinations(n,3)}(index,4) + frac(n)*stroke_source(I,J,K);
                            
                            
                            stroke_source_veins(Destinations(n,1)-numArt-nxnynz) = stroke_source_veins(Destinations(n,1)-numArt-nxnynz) + frac(n)*stroke_source(I,J,K);
                        end
                    end
                end
            end
        end
    end
    stroke_source = stroke_source2;
    stroke_source(stroke_source<error_amt) = 0;
%     sum(sum(sum(stroke_source)))
    stroke_source2 = zeros(size(stroke_source));
end
toc

disp('Adjusting Venous Flowrates (Downstream)')
tic
stroke_source_veins(stroke_source_veins<error_amt) = 0;
stroke_source_veins2 = zeros(size(stroke_source_veins));
while sum(sum(sum(stroke_source_veins))) > 0
    for iv = 1:size(stroke_source_veins,1)
        if stroke_source_veins(iv)
            if ~ismember(iv,OutletPoints)
                Destinations = Tree{iv+nxnynz+numArt};
                if ~isempty(Destinations)
                outflow = sum(Destinations(:,2));
                frac = Destinations(:,2)/outflow;
                for n = 1:size(Destinations,1)

                    if FdotVein(Destinations(n,3),1) <= 0 && FdotVein(Destinations(n,3),2) <= 0
                        FdotVein(Destinations(n,3),:) = FdotVein(Destinations(n,3),:) + frac(n)*stroke_source_veins(iv); % Upstream and Downstream
                    else
                        FdotVein(Destinations(n,3),:) = FdotVein(Destinations(n,3),:) - frac(n)*stroke_source_veins(iv); % Upstream and Downstream
                    end
                    
                    stroke_source_veins2(Destinations(n,1)-numArt-nxnynz) = stroke_source_veins2(Destinations(n,1)-numArt-nxnynz) + frac(n)*stroke_source_veins(iv);
                end
                end
            end
        end
    end
    stroke_source_veins = stroke_source_veins2;
    stroke_source_veins(stroke_source_veins<error_amt) = 0;
%     sum(stroke_source_veins)
    stroke_source_veins2 = zeros(size(stroke_source_veins));
end
toc


disp('Adjusting Arterial Flowrates (Upstream)')
tic
Destinations = Tree2{start_position};

while ~isempty(Destinations)
    Destinations_new = [];
    
    for n = 1:size(Destinations,1)
            
            % Reducing flow down that stream
            if -FdotArt(Destinations(n,3),1) <= 0 && -FdotArt(Destinations(n,3),2) <= 0
                FdotArt(Destinations(n,3),:) = FdotArt(Destinations(n,3),:) + Destinations(n,2); % Upstream and Downstream
            else
                FdotArt(Destinations(n,3),:) = FdotArt(Destinations(n,3),:) - Destinations(n,2); % Upstream and Downstream
            end

            if ~isempty(Tree2{Destinations(n,1)})
                Destination_alter = [Tree2{Destinations(n,1)}(:,1),...
                    (Destinations(n,2)*Tree2{Destinations(n,1)}(:,2)/sum(Tree2{Destinations(n,1)}(:,2))),...
                    Tree2{Destinations(n,1)}(:,3)];
                Destinations_new = [Destinations_new;Destination_alter];
            end
    end
    Destinations = Destinations_new;
%     size(Destinations,1)
end

toc

% redo massFlowrate
measurePerfusion

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
