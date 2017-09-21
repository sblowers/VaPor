function drawVesselsStroke(Vessel)
%drawVessels - Draws a vascular tree using line segments. Black lines indicate unmodified segments and red lines indicate segments modified through vessel generation. (NOTE: Not reccomended if vascular tree contains >20,000 nodes)
%
% Syntax:  drawVessels(Vessel)
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

grid on % Turns grid on.
hold on
view([-45 30])

for Pnt=2:size(Vessel,1); % For all line segments in vessel tree.
    
    PntX=Vessel(Pnt,3); PntY=Vessel(Pnt,4); PntZ=Vessel(Pnt,5);
    % Establish current point locations.
    
    PntConn=Vessel(Pnt,7);
    PntConnX=Vessel(PntConn,3); PntConnY=Vessel(PntConn,4); PntConnZ=Vessel(PntConn,5);
    % Establish connecting point locations.
    
    if (Vessel(Pnt,2)==0), ColourString='r'; else ColourString='k'; end;
    % This colours the vessels according to vessel classification. This
    % defaults to black for vessels from the orginal vessel tree and red
    % for vessels generated or split from vessel generation.
    
    plot3([PntConnY PntY],[PntConnX PntX],[PntConnZ PntZ],ColourString,'LineWidth',2);
    % Plot vessel lines.
    
    drawnow limitrate % Updates figure as it loops.
end
drawnow % Finalises figure.

axis([13 67 5 59 4 58])

set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'ZTickLabel','')

scatter3(Vessel(825,4),Vessel(825,3),Vessel(825,5),50,'ro','LineWidth',2)
view([270,90])