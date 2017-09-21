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


CBF_Full = zeros(size(GM_WM,1),size(GM_WM,2),1);
MTT_calc = zeros(size(GM_WM,1),size(GM_WM,2),1);
CBV_Full = zeros(size(GM_WM,1),size(GM_WM,2),1);

Option_AddVessels = true;

if Option_AddVessels
    addVesselConcentration
end


    if Option_AddVessels
        %     AIF = T_TransientStoreWithVessels(NumDomTot+InletPoints(1),2:end) + T_TransientStoreWithVessels(NumDomTot+InletPoints(2),2:end) + T_TransientStoreWithVessels(NumDomTot+InletPoints(3),2:end);
        AIF = mean([T_TransientStoreWithVessels(NumDomTot+InletPoints(1),2:end); T_TransientStoreWithVessels(NumDomTot+InletPoints(2),2:end); T_TransientStoreWithVessels(NumDomTot+InletPoints(3),2:end)]);
    else
    %     AIF = T_TransientStore(NumDomTot+InletPoints(1),2:end) + T_TransientStore(NumDomTot+InletPoints(2),2:end) + T_TransientStore(NumDomTot+InletPoints(3),2:end);
        AIF = mean([T_TransientStore(NumDomTot+InletPoints(1),2:end); T_TransientStore(NumDomTot+InletPoints(2),2:end); T_TransientStore(NumDomTot+InletPoints(3),2:end)]);
    end
    

% Reducing the timestep from 0.1s to 1s
%     AIF = AIF(linspace(10,300,30));
%     SaveTime = 1;



for X = 1:size(GM_WM,1)
    for Y = 1:size(GM_WM,2)
        for Z = 1:size(GM_WM,3)
        if GM_WM(X,Y,Z)
            
            if Option_AddVessels
                T_Vector = T_TransientStoreWithVessels(GM_WM_Convert(X,Y,Z),2:end) * Porosity2(X,Y,Z);
            else
                T_Vector = T_TransientStore(GM_WM_Convert(X,Y,Z),2:end) * Porosity(X,Y,Z);
            end
            
%             % Reducing the timestep from 0.1s to 1s
%             T_Vector = T_Vector(linspace(10,300,30));
            
%            
            
            
            T_Vector2 = T_Vector;            
            T_Vector2(end+1:2*numel(T_Vector2)-1) = 0;
            Conc = deconv(T_Vector2,AIF).*sum(AIF);  
%             Conc = ifft(fft(T_Vector2)./fft(AIF));
            
            Integral = SaveTime*sum((Conc(1:end-1) + Conc(2:end)))/2;
            MTT = Integral/max(Conc);
            rCBV = (SaveTime*sum((T_Vector(1:end-1) + T_Vector(2:end)))/2)/1;%(SaveTime*sum((AIF(1:end-1) + AIF(2:end)))/2); % Integral of AIF should = 1
            
%             rCBF = max(Conc)/Integral*Porosity2(X,Y,Z)*100*(1000/Rho(X,Y,Z))*60; %  porosity*100*(1000/rho_s) is CBV in ml/100g
%             rCBF = max(Conc)/Integral*2*Porosity(X,Y,Z)*100*(1000/Rho(X,Y,Z))*60; %  porosity*100*(1000/rho_s) is CBV in ml/100g  
            rCBF = rCBV/MTT*100*(1000/Rho(X,Y,Z))*60;
%             rCBF = max(Conc)/Integral*rCBV*100*(1000/mean(Rho(GM_WM)))*60;

%             figure
%             re_conv = conv(C,AIF2)./sum(AIF2);
%             plot(time_vector,Tb_vector,'b',(1:numel(Tb_vector2))*delta_t,Tb_vector2,'g',(1:numel(re_conv))*delta_t,re_conv,'r')
%             legend('C(t)','C(t)_{extended}','C(t)_{reconvoluted}')
            
%             rCBF = GaussNoiseSmooth_insert(Tb_vector,delta_t,20,20,AIF2);
            
            MTT_calc(X,Y,Z) = MTT;
            CBF_Full(X,Y,Z) = rCBF; 
            CBV_Full(X,Y,Z) = rCBV;
            
        else
            CBF_Full(X,Y,Z) = NaN;
            MTT_calc(X,Y,Z) = NaN;
            CBV_Full(X,Y,Z) = NaN;
        end
        end
    end
end

figure
planecut('z',30,CBF_Full)
% planecut('z',10,CBF_Full)
