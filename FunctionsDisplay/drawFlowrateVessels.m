figure
for iv=2:size(Vessel1,1)
    % location of current point
    thisx=Vessel1(iv,3); thisy=Vessel1(iv,4); thisz=Vessel1(iv,5); thisv = abs(Q_art_real(iv,2));
    % find parent location
    parent=Vessel1(iv,7);
    parentx=Vessel1(parent,3); parenty=Vessel1(parent,4); parentz=Vessel1(parent,5); parentv = abs(Q_art_real(iv,1));
    
    [X,Y,Z] = creatingcylinder([parenty,parentx,parentz],[thisy,thisx,thisz],200*Vessel1(parent,6),200*Vessel1(iv,6));
    
    if abs(parentv-thisv) <= 1e-5
        V = ones(1,10)*parentv;
    else
        V = parentv:(thisv-parentv)/9:thisv;
    end
    V = [V;V;V;V;V;V;V;V;V;V;V]';
    
    surf(X,Y,Z,V,'EdgeColor','none','LineStyle','none')
    hold on
end
xlim([1,ny])
ylim([1,nx])
zlim([1,nz])
colorbar
set(gcf,'renderer','zbuffer')