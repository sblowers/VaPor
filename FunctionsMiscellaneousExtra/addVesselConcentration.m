T_TransientStoreWithVessels = T_TransientStore;

Porosity2 = Porosity*VoxelSize^3;

for N = 1:size(T_TransientStoreWithVessels,2)
    T_TransientStoreWithVessels(1:NumGM_WM,N) = T_TransientStoreWithVessels(1:NumGM_WM,N).*((VoxelSize^3)*Porosity(GM_WM));
end

for I = 1:size(GM_WM,1)
    for J = 1:size(GM_WM,2)
        for K = 1:size(GM_WM,3)
            if GM_WM(I,J,K)
                
                if ~isempty(Volume2Vessel1{I,J,K})
                    for M = 1:size(Volume2Vessel1{I,J,K},1)
                        if ~ismember(Volume2Vessel1{I,J,K}(M,1),InletPoints)
                            T_TransientStoreWithVessels(GM_WM_Convert(I,J,K),:) = ...
                                T_TransientStoreWithVessels(GM_WM_Convert(I,J,K),:) + ...
                                (T_TransientStoreWithVessels(NumGM_WM+Volume2Vessel1{I,J,K}(M,1),:)*Vol1(Volume2Vessel1{I,J,K}(M,1))...
                                /size(Vessel2Volume1{Volume2Vessel1{I,J,K}(M,1)},1));
                            
                            Porosity2(I,J,K) = Porosity2(I,J,K) + Vol1(Volume2Vessel1{I,J,K}(M,1))/size(Vessel2Volume1{Volume2Vessel1{I,J,K}(M,1)},1);
                        end
                    end

                end
                
                if ~isempty(Volume2Vessel2{I,J,K})
                    for M = 1:size(Volume2Vessel2{I,J,K},1)
                        T_TransientStoreWithVessels(GM_WM_Convert(I,J,K),:) = ...
                            (T_TransientStoreWithVessels(GM_WM_Convert(I,J,K),:) + ...
                            T_TransientStoreWithVessels(NumGM_WM+NumArt+Volume2Vessel2{I,J,K}(M,1),:)*Vol2(Volume2Vessel2{I,J,K}(M,1))...
                            /size(Vessel2Volume2{Volume2Vessel2{I,J,K}(M,1)},1));
                        
                        Porosity2(I,J,K) = Porosity2(I,J,K) + Vol2(Volume2Vessel2{I,J,K}(M,1))/size(Vessel2Volume2{Volume2Vessel2{I,J,K}(M,1)},1);
                    end
                end
                
            end
        end
    end
end

Porosity2 = Porosity2/VoxelSize^3;
% Porosity2(Porosity2>1)=NaN;

for N = 1:size(T_TransientStoreWithVessels,2)
    T_TransientStoreWithVessels(1:NumGM_WM,N) = T_TransientStoreWithVessels(1:NumGM_WM,N)./((VoxelSize^3)*Porosity2(GM_WM));
end

