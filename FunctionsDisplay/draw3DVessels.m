function draw3DVessels(Vessel,Data,VoxelSize,Limit)
%draw3DVessels - Draws a vascular tree using 3D cylinders. The data is
%coloured based on input Data. A limit is also imposed (default 0.2*VoxelSize) so that
%vessel segments smaller than this diameter are not shown. This can be
%overriden by an input Limit.
%
% Syntax:  draw3DVessels(Vessel,Data,VoxelSize)
%          draw3DVessels(... ,Limit)
%
% Inputs:
%    Vessel - The vessel tree for drawing. Required to be a Nx7 matrix
%    where N is the number of nodes on the vessel tree. Spatial coordinates
%    are in columns 3, 4, & 5 for [x y z] coordinate and column 7 contains
%    the connection node for the vessel tree.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Stephen Blowers - S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
if nargin < 4
    Limit = 0.2*VoxelSize;
end

grid on
colormap jet


for Pnt=2:size(Vessel,1)
    
    if Vessel(Pnt,6) >= Limit
        % location of current point
        PntX=Vessel(Pnt,3); PntY=Vessel(Pnt,4); PntZ=Vessel(Pnt,5);
        % find parent location
        PntConn=Vessel(Pnt,7);
        PntConnX=Vessel(PntConn,3); PntConnY=Vessel(PntConn,4); PntConnZ=Vessel(PntConn,5);
        
        if isempty(Data)
            PntV = 1;
            PntConnV = 1;
        elseif size(Data,2) == 1;
            PntV = Data(Pnt);
            PntConnV = Data(PntConn);
        elseif size(Data,2) == 2;
            PntV = abs(Data(Pnt,2));
            PntConnV = abs(Data(Pnt,1));
        end
        
        [X,Y,Z] = createCylinder([PntConnY,PntConnX,PntConnZ],[PntY,PntX,PntZ],Vessel(PntConn,6)/VoxelSize,Vessel(Pnt,6)/VoxelSize);
        
        if abs(PntConnV-PntV) <= 1e-5
            V = ones(1,10)*PntConnV;
        else
            V = PntConnV:(PntV-PntConnV)/9:PntV;
        end
        V = [V;V;V;V;V;V;V;V;V;V;V]';
        
        surf(X,Y,Z,V,'EdgeColor','none','LineStyle','none')
        drawnow limitrate
        hold on
    end
end
drawnow

view([-30 30])
light('Position',[0,1,1])
material shiny

if ~isempty(Data)
    colorbar
end
set(gcf,'renderer','zbuffer')