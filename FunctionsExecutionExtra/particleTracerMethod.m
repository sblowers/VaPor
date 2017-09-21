disp('%%%%%%% Running Particle Tracer Method %%%%%%%%%%%%')
disp('Creating Network Trees')
tic

%%%%%% Setting Up Reference Lists for DomTot %%%%
NumArt = size(Vessel1,1);
NumVein = size(Vessel2,1);
NumGM_WM = numel(GM_WM(GM_WM));
NumDomTot = NumGM_WM;
RowConvert = find(GM_WM);  % goes from row number to I,J,K
GM_WM_Convert = zeros(size(GM_WM));
GM_WM_Convert(RowConvert) =  (1:numel(GM_WM(GM_WM))); % goes from I,J,K to row number
% The sparse matrix to solve only needs those voxels that are within
% DomTot. These lists help convert sparse matrix row number to location in
% DomTot (given by I,J,K) and from locations in DomTot to row numbers in
% the sparse matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MTT = (Porosity*VoxelSize^3)./(MassFlow/Rho_b); % s
% MTT = 1./(MassFlow/rho./(porosity*h^3))*3; % s
% Guess for mean transit time for voxels given porosity.

Tree = cell(NumGM_WM+NumArt+NumVein,1);
% data stored in [location, probability to travel to that location (based on flow splits),time to travel]


% find downstream flows and weights for arteries
for ArtLoop = 1:NumArt
    
    if SplitArt(ArtLoop,1)~=0
        for N = 1:size(Vessel2Volume1{ArtLoop},1)
                if Vessel2Volume1{ArtLoop}(N,4) ~= 0;
                    Tree{ArtLoop+NumGM_WM}(end+1,:) = [GM_WM_Convert(Vessel2Volume1{ArtLoop}(N,1),Vessel2Volume1{ArtLoop}(N,2),Vessel2Volume1{ArtLoop}(N,3))...
                        ,SplitArt(ArtLoop,2)*Vessel2Volume1{ArtLoop}(N,4),0];
                    Tree{Vessel1(ArtLoop,7)+NumGM_WM}(end+1,:) = [GM_WM_Convert(Vessel2Volume1{ArtLoop}(N,1),Vessel2Volume1{ArtLoop}(N,2),Vessel2Volume1{ArtLoop}(N,3))...
                        ,(1-SplitArt(ArtLoop,2))*Vessel2Volume1{ArtLoop}(N,4),0];
                end
        end
        
        
    else
    
    
    if FdotArt(ArtLoop,1) <= 0 
        % probablility of travelling down artery length
        if abs(FdotArt(ArtLoop,1)) > 0
            Tree{ArtLoop+NumGM_WM}(end+1,:) = [Vessel1(ArtLoop,7)+NumGM_WM,abs(FdotArt(ArtLoop,1)),-Vol1(ArtLoop)/(abs(FdotArt(ArtLoop,1))/Rho_b)];
        end
        
        % probablility of transferring domain on artery length
        if ~isempty(Vessel2Volume1{ArtLoop})
            for N = 1:size(Vessel2Volume1{ArtLoop},1)
                if Vessel2Volume1{ArtLoop}(N,4) ~= 0;
                    Tree{ArtLoop+NumGM_WM}(end+1,:) = [GM_WM_Convert(Vessel2Volume1{ArtLoop}(N,1),Vessel2Volume1{ArtLoop}(N,2),Vessel2Volume1{ArtLoop}(N,3))...
                        ,Vessel2Volume1{ArtLoop}(N,4),0];
                end
            end
        end
        
    elseif FdotArt(ArtLoop,2) >= 0 
        
        % probablility of travelling down artery length
        if abs(FdotArt(ArtLoop,2)) > 0 
            Tree{Vessel1(ArtLoop,7)+NumGM_WM}(end+1,:) = [ArtLoop+NumGM_WM,abs(FdotArt(ArtLoop,2)),Vol1(ArtLoop)/(abs(FdotArt(ArtLoop,2))/Rho_b)];
        end
        
        % probablility of transferring domain on artery length
        if ~isempty(Vessel2Volume1{ArtLoop})
            for N = 1:size(Vessel2Volume1{ArtLoop},1)
                if Vessel2Volume1{ArtLoop}(N,4) ~= 0;
                    Tree{Vessel1(ArtLoop,7)+NumGM_WM}(end+1,:) = [GM_WM_Convert(Vessel2Volume1{ArtLoop}(N,1),Vessel2Volume1{ArtLoop}(N,2),Vessel2Volume1{ArtLoop}(N,3))...
                        ,Vessel2Volume1{ArtLoop}(N,4),0];
                end
            end
        end
        
    end
    
    end
    
end

% find downstream flows and weights for domain
for DomainLoop = 1:NumGM_WM
    DownTemp = [];
    
    [I,J,K] = ind2sub(size(GM_WM),RowConvert(DomainLoop));

    Ax = VoxelSize^2; Ay = VoxelSize^2; Az = VoxelSize^2;
    if U(I,J,K)<0, DownTemp(end+1,:) = [GM_WM_Convert(I-1,J,K),Ax*Rho_b*abs(U(I,J,K)),-MTT(I,J,K)]; end % MTT for domain already calculated in massFlowrate
    if U(I+1,J,K)>0, DownTemp(end+1,:) = [GM_WM_Convert(I+1,J,K),Ax*Rho_b*abs(U(I+1,J,K)),-MTT(I,J,K)]; end
    if V(I,J,K)<0, DownTemp(end+1,:) = [GM_WM_Convert(I,J-1,K),Ay*Rho_b*abs(V(I,J,K)),-MTT(I,J,K)]; end
    if V(I,J+1,K)>0, DownTemp(end+1,:) = [GM_WM_Convert(I,J+1,K),Ay*Rho_b*abs(V(I,J+1,K)),-MTT(I,J,K)]; end
    if W(I,J,K)<0, DownTemp(end+1,:) = [GM_WM_Convert(I,J,K-1),Az*Rho_b*abs(W(I,J,K)),-MTT(I,J,K)]; end
    if W(I,J,K+1)>0, DownTemp(end+1,:) = [GM_WM_Convert(I,J,K+1),Az*Rho_b*abs(W(I,J,K+1)),-MTT(I,J,K)]; end

    if ~isempty(Volume2Vessel2{I,J,K})
        for N = 1:size(Volume2Vessel2{I,J,K},1)
            DownNode = Volume2Vessel2{I,J,K}(N,1);
            if SplitVein(DownNode,1)~=0
                DownTemp(end+1,:) = [DownNode+NumArt+NumGM_WM,SplitVein(N,2)*abs(Volume2Vessel2{I,J,K}(N,2)),-MTT(I,J,K)];
                DownTemp(end+1,:) = [Vessel2(DownNode,7)+NumArt+NumGM_WM,(1-SplitVein(N,2))*abs(Volume2Vessel2{I,J,K}(N,2)),-MTT(I,J,K)];
            else
                if FdotVein(DownNode,2) <= 0
                    DownNode = Vessel2(DownNode,7);
                end
                DownTemp(end+1,:) = [DownNode+NumArt+NumGM_WM,abs(Volume2Vessel2{I,J,K}(N,2)),-MTT(I,J,K)];
                
            end
        end
    end
    
    Tree{DomainLoop} = DownTemp;
    
end

% find downstream flows and weights for veins
for VeinLoop = 1:NumVein  
    
    if FdotVein(VeinLoop,2) < 0 % && abs(FdotVein(VeinLoop,2)) > ErrorAmt
        Tree{VeinLoop+NumArt+NumGM_WM}(end+1,:) = [Vessel2(VeinLoop,7)+NumArt+NumGM_WM,abs(FdotVein(VeinLoop,2)),-Vol2(VeinLoop)/(abs(FdotVein(VeinLoop,2))/Rho_b)];
    elseif FdotVein(VeinLoop,1) > 0 % && abs(FdotVein(VeinLoop,1)) > ErrorAmt
        Tree{Vessel2(VeinLoop,7)+NumArt+NumGM_WM}(end+1,:) = [VeinLoop+NumArt+NumGM_WM,abs(FdotVein(VeinLoop,1)),Vol2(VeinLoop)/(abs(FdotVein(VeinLoop,1))/Rho_b)];
    end
    
end



for N = 1:size(Tree,1)
    if ~isempty(Tree{N})
        Tree{N}(:,2) = Tree{N}(:,2)/sum(Tree{N}(:,2));
    end
end
toc


% travelLog gives position and then time spent in that position
% figure
% hold on
% grid on

% timelist = zeros(num_points_generated,4);

% initialise concentration measurements
Timestep = 0.1;
SaveTime = Timestep;
TotalTime = 30;
NoTimesteps = ceil(TotalTime/Timestep);
TotalTime = NoTimesteps*Timestep;

InjectionLength = 1;
T_TransientStore = zeros(size(Tree,1),TotalTime/SaveTime);

disp('Simulating Particles')
tic


% parfor par_loop = 1:8
% for par_loop = 1:4
%     con_domains_temp = zeros(size(Tree,1),max_time/delta_t + 1);

Option_Streamlines = false;
StreamlineFigureTrigger = true;

timelist = zeros(ParticlesGenerated,4);

for npnts = 1:ParticlesGenerated
    if rem(npnts,10000) == 0
        disp([num2str(npnts) ' Iterations Completed'])
    end
    
    AllocatedRows = 1000;
    travelLog = zeros(AllocatedRows,2);
    rowcount=1;
    %             tic
    % point path function
    % choose inlet
    location = InletPoints(find(cumsum(InletFlows/sum(InletFlows))>=rand,1,'first'))+NumGM_WM;
    travelLog(rowcount,:) = [location(1),0];
    rowcount=rowcount+1;
    % update location
    
    
    while  ~ismember(location,OutletPoints+NumArt+NumGM_WM)
        location_matrix = Tree{location};
        location_new = location_matrix(find(cumsum(location_matrix(:,2))>=rand,1,'first'),[1,3]);
        
        %             sum_col4 = sum(location_matrix(:,4));
        %             if sum_col4==0, sum_col4=1; end
        %             [max_val,max_ind] = max(location_matrix(:,2)-location_matrix(:,4)/sum_col4);
        %             location_new = location_matrix(max_ind,[1,3]);
        %             Tree{location}(max_ind,4) = Tree{location}(max_ind,4) + 1;
        
        travelLog(rowcount,:) = [location_new(1,1),0];
        
        % update time
                    if location_new(1,1) <= NumGM_WM
%                         fudged_MTT = location_new(1,2) + ((2*rand)-1)*location_new(1,2);
                        fudged_MTT = normrnd(location_new(1,2),abs(location_new(1,2)/4));
                    else
                        fudged_MTT = location_new(1,2);
                    end
%         fudged_MTT = location_new(1,2);
        if  location_new(2) < 0
            travelLog(rowcount-1,2) = travelLog(rowcount-1,2) + abs(fudged_MTT);
        else
            travelLog(rowcount,2) = abs(fudged_MTT);
        end
        
        %             travelLog(rowcount,2) = travelLog(rowcount,2) + ((2*rand)-1)*travelLog(rowcount,2);
        
        location = location_new(1,1);
        rowcount=rowcount+1;
    end
    
    if rowcount > AllocatedRows
        disp(['TravelLog larger than allocated rows ( ' num2str(rowcount) ' > ' num2str(AllocatedRows) ')'])
    end
    
    travelLog = travelLog(logical(travelLog(:,1)),:);
    %         toc
    
    % draw streamline
    if Option_Streamlines && npnts <= 1000
        if StreamlineFigureTrigger 
            figure; hold on; grid on; axis([0 72 0 72 0 72]); colorbar; view([-30,30]);
            StreamlineFigureTrigger = false;
        end
    
            loc = zeros(size(travelLog,1),3);
        for n = 1:size(travelLog,1)
            
            
            
            if travelLog(n,1) <=NumGM_WM 
                [loc(n,1),loc(n,2),loc(n,3)] = ind2sub(size(GM_WM),RowConvert(travelLog(n,1)));
                loc(n,:) = loc(n,:) + [rand-0.5,rand-0.5,rand-0.5];
            elseif travelLog(n,1)<=NumGM_WM+NumArt
                loc(n,:) = Vessel1(travelLog(n,1)-NumGM_WM,3:5);
                loc(n,:) = loc(n,:)+ sph2cart(asin(2*rand-1), 2*pi*rand, (rand-0.5)*Vessel1(travelLog(n,1)-NumGM_WM,6)/VoxelSize);
            else
                loc(n,:) = Vessel2(travelLog(n,1)-NumGM_WM-NumArt,3:5);
                loc(n,:) = loc(n,:)+ sph2cart(asin(2*rand-1), 2*pi*rand, (rand-0.5)*Vessel2(travelLog(n,1)-NumGM_WM-NumArt,6)/VoxelSize);
            end
        end
        
        loc = [loc(:,2),loc(:,1),loc(:,3)]; % Switch x and y because of how matlab plots
        
        ArtEnd = find(travelLog(:,1)<NumGM_WM,1,'first')-1;
        GM_WM_End = find(travelLog(:,1)<NumGM_WM,1,'last')+1;
        

%         plot3(loc(1:ArtEnd,1),loc(1:ArtEnd,2),loc(1:ArtEnd,3),'r','LineWidth',1)
%         plot3(loc(GM_WM_End:end,1),loc(GM_WM_End:end,2),loc(GM_WM_End:end,3),'b','LineWidth',1)
%         fnplt(cscvn(loc(ArtEnd:GM_WM_End,:)'),'k',1)

        % Plot coloured by MTT with surfaces
        MTT_Total = cumsum(travelLog(:,2));
               
%         SplinePoints = fnplt(cscvn(loc(ArtEnd:GM_WM_End,:)'));
%         MTTSplinePoints = interp1(1:numel(ArtEnd:GM_WM_End),MTT_Total(ArtEnd:GM_WM_End),linspace(1,numel(ArtEnd:GM_WM_End),size(SplinePoints,2)));
        SplinePoints = fnplt(cscvn(loc'));
        MTTSplinePoints = interp1(1:numel(MTT_Total),MTT_Total,linspace(1,numel(MTT_Total),size(SplinePoints,2)));
%         

%         surface([loc(1:ArtEnd,1),loc(1:ArtEnd,1)]',[loc(1:ArtEnd,2),loc(1:ArtEnd,2)]',...
%             [loc(1:ArtEnd,3),loc(1:ArtEnd,3)]',[MTT_Total(1:ArtEnd) MTT_Total(1:ArtEnd)]',...
%             'facecol','no', 'edgecol','interp','linew',1);
%         surface([loc(GM_WM_End:end,1),loc(GM_WM_End:end,1)]',[loc(GM_WM_End:end,2),loc(GM_WM_End:end,2)]',...
%             [loc(GM_WM_End:end,3),loc(GM_WM_End:end,3)]',[MTT_Total(GM_WM_End:end),MTT_Total(GM_WM_End:end)]',...
%             'facecol','no', 'edgecol','interp','linew',1);

        surface([SplinePoints(1,:);SplinePoints(1,:)],[SplinePoints(2,:);SplinePoints(2,:)],...
            [SplinePoints(3,:);SplinePoints(3,:)],[MTTSplinePoints;MTTSplinePoints],...
            'facecol','no', 'edgecol','interp','linew',1);
        colorbar
        view([-30,30])
    end
    
    
    timelist(npnts,1) = sum(travelLog(:,2));
    timelist(npnts,2) = sum(travelLog(1:find(travelLog(:,1)<NumGM_WM,1,'first')-1,2));
    timelist(npnts,3) = sum(travelLog(find(travelLog(:,1)<NumGM_WM,1,'first'):find(travelLog(:,1)<=NumGM_WM,1,'last'),2));
    timelist(npnts,4) = sum(travelLog(find(travelLog(:,1)<=(NumGM_WM),1,'last')+1:end,2));
    
    %
    % tic
    cumm_time = cumsum(travelLog(:,2));
    %         con_domains(travelLog(1,1),1) = con_domains(travelLog(1,1),1) + 1;
    
    for N = 1:(TotalTime/Timestep)
        t = N*Timestep;
        
        if t-InjectionLength <= cumm_time(end);
            
            %             bar1 = find(cumm_time>=t,1,'first');
            %             con_domains(travelLog(bar1,1),n+1) = con_domains(travelLog(bar1,1),n+1) + 1;
            
            bar1 = find(cumm_time>=t,1,'first'); if isempty(bar1), bar1 = numel(cumm_time); end
            bar2 = find(cumm_time>=(t-InjectionLength),1,'first');
            %             con_domains(travelLog(bar1:bar2,1),n+1) = con_domains(travelLog(bar1:bar2,1),n+1) + 1;
            
            if bar1 == bar2
                if bar1 == numel(cumm_time)
                    T_TransientStore(travelLog(bar1,1),N) = T_TransientStore(travelLog(bar1,1),N) + min((cumm_time(end)-(t-InjectionLength)),InjectionLength);
                else
                    T_TransientStore(travelLog(bar1,1),N) = T_TransientStore(travelLog(bar1,1),N) + min(t,InjectionLength);
                end
            else
                T_TransientStore(travelLog(bar1,1),N) = T_TransientStore(travelLog(bar1,1),N) + t-cumm_time(bar1-1);
                T_TransientStore(travelLog(bar2,1),N) = T_TransientStore(travelLog(bar2,1),N) + cumm_time(bar2)-max((t-InjectionLength),0);
                T_TransientStore(travelLog(bar2+1:bar1-1,1),N) = T_TransientStore(travelLog(bar2+1:bar1-1,1),N) + cumm_time(bar2+1:bar1-1)-cumm_time(bar2:bar1-2);
            end
            
        end
    end
    
    %         cumm_time = cumsum(travelLog(:,2));
    %         con_domains_temp(travelLog(1,1),1) = con_domains_temp(travelLog(1,1),1) + 1;
    %
    %         for n = 1:(max_time/delta_t)
    %             t = n*delta_t;
    %
    %             if t-inj_length <= cumm_time(end);
    %
    %                 %             bar1 = find(cumm_time>=t,1,'first');
    %                 %             con_domains(travelLog(bar1,1),n+1) = con_domains(travelLog(bar1,1),n+1) + 1;
    %
    %                 bar1 = find(cumm_time>=t,1,'first'); if isempty(bar1), bar1 = numel(cumm_time); end
    %                 bar2 = find(cumm_time>=(t-inj_length),1,'first');
    %                 %             con_domains(travelLog(bar1:bar2,1),n+1) = con_domains(travelLog(bar1:bar2,1),n+1) + 1;
    %
    %                 if bar1 == bar2
    %                     if bar1 == numel(cumm_time)
    %                         con_domains_temp(travelLog(bar1,1),n+1) = con_domains_temp(travelLog(bar1,1),n+1) + min((cumm_time(end)-t),inj_length);
    %                     else
    %                         con_domains_temp(travelLog(bar1,1),n+1) = con_domains_temp(travelLog(bar1,1),n+1) + min(t,inj_length);
    %                     end
    %                 else
    %                     con_domains_temp(travelLog(bar1,1),n+1) = con_domains_temp(travelLog(bar1,1),n+1) + t-cumm_time(bar1-1);
    %                     con_domains_temp(travelLog(bar2,1),n+1) = con_domains_temp(travelLog(bar2,1),n+1) + cumm_time(bar2)-max((t-inj_length),0);
    %                     con_domains_temp(travelLog(bar2+1:bar1-1,1),n+1) = con_domains_temp(travelLog(bar2+1:bar1-1,1),n+1) + cumm_time(bar2+1:bar1-1)-cumm_time(bar2:bar1-2);
    %                 end
    %
    %             end
    %         end
    % toc
end
%
%     con_domains = con_domains + con_domains_temp;
% end

% T_TransientStore = [zeros(size(Tree,1),1), [T_TransientStore(NumArt+1:NumArt+NumGM_WM,:);T_TransientStore(1:NumArt,:);T_TransientStore(NumArt+NumGM_WM+1:end,:)]];
T_TransientStore_Save = T_TransientStore;

T_TransientStore = T_TransientStore/ParticlesGenerated;
% 
for N = 1:size(T_TransientStore,2)
    T_TransientStore(1:NumGM_WM,N) = T_TransientStore(1:NumGM_WM,N)./((VoxelSize^3)*Porosity(GM_WM));
    T_TransientStore(NumGM_WM+2:NumGM_WM+NumArt,N) = T_TransientStore(NumGM_WM+2:NumGM_WM+NumArt,N)./Vol1(2:end);
    T_TransientStore(NumGM_WM+1,N) = T_TransientStore(NumGM_WM+2,N);
    T_TransientStore(NumGM_WM+NumArt+2:end,N) = T_TransientStore(NumGM_WM+NumArt+2:end,N)./Vol2(2:end);
    T_TransientStore(NumGM_WM+NumArt+1,N) = T_TransientStore(NumGM_WM+NumArt+2,N);
end

T_TransientStore = [zeros(size(Tree,1),1), T_TransientStore];
% T_TransientStore(NumGM_WM+InletPoints,:) = [ones(1,round(InjectionLength/Timestep)+1), zeros(1,round(TotalTime/Timestep-(InjectionLength/Timestep)));...
%     ones(1,round(InjectionLength/Timestep)+1), zeros(1,round(TotalTime/Timestep-(InjectionLength/Timestep)));...
%     ones(1,round(InjectionLength/Timestep)+1), zeros(1,round(TotalTime/Timestep-(InjectionLength/Timestep)))];

AIF = mean([T_TransientStore(NumDomTot+InletPoints(1),2:end); T_TransientStore(NumDomTot+InletPoints(2),2:end); T_TransientStore(NumDomTot+InletPoints(3),2:end)]);
T_TransientStore = T_TransientStore/max(AIF);


toc

% addVesselConcentration

plotPerfusionTracer