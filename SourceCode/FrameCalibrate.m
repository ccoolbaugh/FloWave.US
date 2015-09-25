function Conversion = FrameCalibrate(I)
%FRAMECALIBRATE Pixel conversion to units of interest
%   
%   Conversion = FRAMECALIBRATE(I) gets image of the scale axis (I) and
%   returns the pixel scaling factor. The user selects 2 positions on the
%   image scale using the function GIPNUTC to define a range for a scale
%   factor. The scale factor is entered by the user as an input and the
%   number of pixels are converted to a unit of interest. 
%
%   The positions selected on the scale should be at extremes to improve
%   the accuracy of the scaling factor. 
%

% Crystal Coolbaugh
% April 16, 2015
% Copyright 2015 Crystal Coolbaugh


% Display Scale
figure;
image(I);

title('Scale Calibration', 'FontSize', 16)

disp('Select the scale extremes (e.g. 0 cm/s and 80 cm/s)')
[Xpos, Ypos] = ginputc(2, 'Color','r','LineWidth',3);

Direction = menu('Choose the scale direction','Vertical','Horizontal');

if Direction == 1
    PixDiff = abs(diff(Ypos));
else
    PixDiff = abs(diff(Xpos));
end

ScaleFactor = input('Enter the scale conversion factor (e.g. 10 cm/s = 10): ');
Conversion = ScaleFactor/PixDiff;

clear Xpos Ypos
close;
