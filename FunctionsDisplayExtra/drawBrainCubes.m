figure
hold on
Color1 = [240,200,201]./255;
% Color2 = Color1;
Color2 = [0,0,1];

[I,J,K] = ind2sub(size(Borders),find(Borders==1));
CornersX = zeros(4,numel(I));
CornersY = zeros(4,numel(J));
CornersZ = zeros(4,numel(K));
Color = zeros(1,numel(I),3);
for M = 1:6
    for N = 1:numel(I)
        switch M
            case 1
                CornersX(:,N) = [I(N)-0.5;I(N)-0.5;I(N)-0.5;I(N)-0.5];
                CornersY(:,N) = [J(N)+0.5;J(N)+0.5;J(N)-0.5;J(N)-0.5];
                CornersZ(:,N) = [K(N)+0.5;K(N)-0.5;K(N)-0.5;K(N)+0.5];
            case 2
                CornersX(:,N) = [I(N)+0.5;I(N)+0.5;I(N)+0.5;I(N)+0.5];
                CornersY(:,N) = [J(N)+0.5;J(N)+0.5;J(N)-0.5;J(N)-0.5];
                CornersZ(:,N) = [K(N)+0.5;K(N)-0.5;K(N)-0.5;K(N)+0.5];
            case 3
                CornersX(:,N) = [I(N)+0.5;I(N)+0.5;I(N)-0.5;I(N)-0.5];
                CornersY(:,N) = [J(N)-0.5;J(N)-0.5;J(N)-0.5;J(N)-0.5];
                CornersZ(:,N) = [K(N)+0.5;K(N)-0.5;K(N)-0.5;K(N)+0.5];
            case 4
                CornersX(:,N) = [I(N)+0.5;I(N)+0.5;I(N)-0.5;I(N)-0.5];
                CornersY(:,N) = [J(N)+0.5;J(N)+0.5;J(N)+0.5;J(N)+0.5];
                CornersZ(:,N) = [K(N)+0.5;K(N)-0.5;K(N)-0.5;K(N)+0.5];
            case 5
                CornersX(:,N) = [I(N)+0.5;I(N)+0.5;I(N)-0.5;I(N)-0.5];
                CornersY(:,N) = [J(N)+0.5;J(N)-0.5;J(N)-0.5;J(N)+0.5];
                CornersZ(:,N) = [K(N)-0.5;K(N)-0.5;K(N)-0.5;K(N)-0.5];
            case 6
                CornersX(:,N) = [I(N)+0.5;I(N)+0.5;I(N)-0.5;I(N)-0.5];
                CornersY(:,N) = [J(N)+0.5;J(N)-0.5;J(N)-0.5;J(N)+0.5];
                CornersZ(:,N) = [K(N)+0.5;K(N)+0.5;K(N)+0.5;K(N)+0.5];
        end
        
        if K(N) > 42
            Color(1,N,:) = Color2;
        else
            Color(1,N,:) = Color1;
        end
    end
    patch(CornersX,CornersY,CornersZ,Color)
    
end

axis([0 72 1 73 0 72])


plot3(ones(40,1)*30,ones(40,1)*36,linspace(30,72,40),'k','LineWidth',4)
plot3(ones(40,1)*30,ones(40,1)*30,linspace(30,72,40),'k','LineWidth',4)
plot3(ones(40,1)*30,ones(40,1)*24,linspace(30,72,40),'k','LineWidth',4)

grid on

set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'ZTickLabel','')

view([100 25])