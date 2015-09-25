function StartFrame = FigureLoopingFcn(VidObj,FramesPerPeriod)
%FIGURELOOPINGFCN Anitmates figure to determine video start frame.
%   FIGURELOOPINGFCN creates an animated figure of individual video frames.
%
%   StartFrame = FIGURELOOPINGFCN(VidObj, FramesPerPeriod) gets the video
%   object and the number of frames per pulse wave update and returns the
%   frame number for the start of the video as indicated by the user. The
%   user stops the animated frame scroll by pressing a key (callback to
%   myKeyPressFcn). The frame identified by the user is saved as the
%   initial start frame for the video. 
%

% Crystal Coolbaugh
% July 28, 2015
% Copyright 2015 Crystal Coolbaugh

global ScrollStop

ScrollStop = 0; %Flag to stop video frame loop
k = FramesPerPeriod;  %Frame count

fig_start = figure;
set(fig_start,'KeyPressFcn',@myKeyPressFcn)

while ~ScrollStop
    vidframe = read(VidObj,k);
    image(vidframe);
    title(['Frame Number: ' int2str(k)]);
    pause(0.2);
    
    k = k+1;
end

StartFrame = k - 1;

function myKeyPressFcn(hObject, event)

global ScrollStop

ScrollStop = 1;




 
