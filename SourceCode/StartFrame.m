function [StartFrame] = StartFrame(USObj,FramesPerPeriod)
%%% Identify Video Start Frame
% The ultrasound video should be edited with an external program so that
% the pulse wave data on the initial start frame of the video represents a
% complete update of data. To account for errors with editing, the user is
% presented an animation of the start of the video with the individual
% frame numbers. The user can press the  spacebar to stop the animation at
% the correct start frame. 
 
disp('An animation of the video frames will be displayed.')
disp('Press the spacebar to stop the video when the update of the pulse wave data is at the far right side of the screen (i.e. a full pulse wave update).')
disp('Paused. Press any key to continue.')
pause;

StartCheck = 0;

while ~StartCheck
    StartFrame = FigureLoopingFcn(USObj,FramesPerPeriod);
    
    % Check Frame Selection
    disp('Start Frame Selected: ')
    StartFrame
    
    StartTest = input('Was the start frame selected correctly (Y/N)?', 's');
    if lower(StartTest) == 'y'
        disp('Advancing to TAMean Extraction')
        close all;
        StartCheck = 1;
    end
end


end

