function planecutTracer(Dir,Val,Data)
%planecut - Takes a slice of an input Data at specificed location using the inbuilt Matlab 'slice' function. Removes any lines segmenting the data and orientates for better viewing.
%
% Syntax:  planecut(Dir,Val,Data)
%
% Inputs:
%    Dir - Direction for data to be sliced. A string value that can either
%    be 'x', 'y', or 'z'.
%    Val - Location within data to take a slice. Can be a number between 1
%    and size of Data in direction Dir. Does not necessarily have to be an
%    integer value.
%    Data - 3D matrix to be sliced.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Stephen Blowers - S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017


if numel(size(Data))==2 || (numel(size(Data))==3 && any(size(Data)==1))
    Data = squeeze(Data);
    H = surf(ones(size(Data)),Data); % If data is 2D.
    view([270 90]) % Sets view to see the sliced data.
else
    
    if strcmp(Dir,'x')
        H = slice(Data,Val,[],[]); % Takes a slice of the data in the X direction.
        view([-90 0]) % Sets view to see the sliced data.
    end
    
    if strcmp(Dir,'y')
        H = slice(Data,[],Val,[]); % Takes a slice of the data in the X direction.
        view([0 0]) % Sets view to see the sliced data.
    end
    
    if strcmp(Dir,'z')
        H = slice(Data,[],[],Val); % Takes a slice of the data in the X direction.
        view([270 90]) % Sets view to see the sliced data.
    end
end

H.EdgeColor = 'none';
H.LineStyle = 'none';
% Turns off any edges that appear within the slice of data.

alpha(H,'color')

AxisLength = max(size(Data));
axis([0 AxisLength 0 AxisLength 0 AxisLength])
% Sets the axis to be square.

colorbar % Shows colourbar.
colormap([ones(100,1) linspace(1,0,100)' linspace(1,0,100)'])
caxis([0 1])
end