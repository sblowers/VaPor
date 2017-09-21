VideoFilename = 'VideoTest.mp4';
VideoLength = 15; % Length of video [s]
NoFrames = floor(VideoLength/SaveTime + 1); % Number of frames in video [].


WriterObj = VideoWriter(VideoFilename,'MPEG-4');
WriterObj.FrameRate = 1/SaveTime;
% writerObj.Quality = 100;

open(WriterObj);

if exist('F','var')
    clearvars('F')
end
F(NoFrames) = struct('cdata',[],'colormap',[]);

for Frame = 1:NoFrames
    tic
    disp(['Drawing Frame: ' num2str(Frame) '/' num2str(NoFrames)])
    
    % Draw Arteries
    T_Art = T_TransientStore(NumDomRows+1:NumDomRows+VesselRow,Frame);
    
    Fig = figure;
    set(Fig,'units','normalized','outerposition',[0 0 1 1],'Visible','off');
%     set(Fig,'units','normalized','outerposition',[0 0 1 1]);
%     draw3DVesselsTracer(Vessel1,T_Art,VoxelSize,0.00001)
    draw3DVesselsTracer(Vessel1,T_Art,VoxelSize)
%     F(Frame) = getframe(Fig);
%     close(Fig)
%         
%     % Draw Veins
    T_Vein = T_TransientStore(NumDomRows+VesselRow+1:end,Frame);
%     
%     Fig = figure;
%     set(Fig,'units','normalized','outerposition',[0 0 1 1],'Visible','off');
%     draw3DVesselsTracer(Vessel2,T_Vein,VoxelSize,0.00001)
    draw3DVesselsTracer(Vessel2,T_Vein,VoxelSize)
    
    caxis([0 1])
    set(gcf,'renderer','zbuffer')
    
    F(Frame) = getframe(Fig);
    close(Fig)
        
%     Draw Planecut
%     Tt = zeros(size(DomTot));
%     Tt(RowConvert) = T_TransientStore(1:NumDomRows,Frame);
%     
%     Fig = figure;
%     set(Fig,'units','normalized','outerposition',[0 0 1 1],'Visible','off');
%     planecutTracer('z',30,Tt)
%     F(Frame) = getframe(Fig);
%     close(Fig)
    toc
end

writeVideo(WriterObj,F);
close(WriterObj)