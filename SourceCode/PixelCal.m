function [VelCon,TimeCon,DistCon] = PixelCal(OneFrame,PlatSetFile)
%%% Pixel Calibration (FrameCalibrate Function)
% Determine the factor needed to convert pixels to the quantity of
% interest: velocity, time, and distance. 
% Position of the scales may need to be adjusted for different US
% equipment. Inspect a single image frame and adjust the zoomed portion.

load(PlatSetFile)

ScaleCheck = 0;
disp('Use the cursor to select the scale range (lower & upper extremes).')

while ScaleCheck == 0
    % Set Velocity Conversion Factor - Zoom to Velocity Scale
    VelScale = OneFrame(VelY(1):VelY(2), VelX(1):VelX(2), :);
    VelCon = FrameCalibrate(VelScale);
    
    % Set Time Conversion Factor - Zoom to Time Scale
    TimeScale = OneFrame(TimeY(1):TimeY(2),TimeX(1):TimeX(2), :);
    TimeCon = FrameCalibrate(TimeScale);
    
    % Set Distance Conversion Factor - Zoom to Depth Scale
    DistScale = OneFrame(DistY(1):DistY(2),DistX(1):DistX(2), :);
    DistCon = FrameCalibrate(DistScale);
    
    % Calibration Check
    CalCheck = input('Were calibrations performed correctly?(Y/N)? ','s');
    
    if lower(CalCheck) == 'y'
        ScaleCheck = 1;
    end
end


end

