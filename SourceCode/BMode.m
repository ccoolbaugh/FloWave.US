%% BMode
% Purpose: This code was written to analyze ultrasound (US) B-mode images
% from a digital video recording. The code uses an automated diameter 
% subroutine to determine a blood vessel diameter. 
%
% Inputs: Ultrasound Video (type - '.avi' or '.mov')
% Outputs: Vessel diameter (cm)
%
% Functions: FrameCalibrate, VesselROI, AutoDiameter, imgaussian (written
% by Dirk-Jan Kroon), ginputc (written by Jiro Doke) - See associated
% license agreements for copyright information.
%
% Requires MATLAB R2012 or newer.
%
% Crystal Coolbaugh
% August 4, 2015
% Copyright 2015 Crystal Coolbaugh

%% Format Workspace
clear
format compact
clc
close all
disp('PLEASE FOLLOW THE ON-SCREEN PROMPTS TO PROCESS THE ULTRASOUND DATA.')

%% Identify Folder with Video
% Add video folder to the current path
disp('Choose the file directory that contains the US video.');
DirName=uigetdir;   %Open browser to identify the directory path
addpath(DirName);
ls(DirName); %List file contents for easy copy-past of filename

%% Import US Video Data

% Enter the Video Filename - Can Copy/Paste from List
VideoName = input('Type the video filename and extension (e.g. USVideo.avi): ','s');

%Construct an object using the VideoReader Function
USObj = VideoReader(VideoName);

% ACQUIRE VIDEO SETTINGS: May Vary for Different US or DVD equipment
% Set Image Frame Size 
VidHeight = USObj.Height;
VidWidth = USObj.Width;

% Total Number of Image Frames
nFrames = USObj.NumberofFrames;

% Video Sampling Rate
VidSamRate = USObj.FrameRate;

% Video End Time
Duration = USObj.Duration;

%% Pixel Calibration (FrameCalibrate Function)
% Determine the factor needed to convert pixels to the quantity of
% interest: distance. Position of the scales may need to be adjusted for
% different US equipment or screen resolutions. 

%Identify position of the distance scale
OneFrame = read(USObj, 1);
image(OneFrame); colormap gray
title('Define a large ROI around the distance scale (click 1: upper left corner; click 2: lower right corner)');
[X,Y] = ginputc(2, 'Color', 'r', 'LineWidth', 2);
ScaleX = floor(X);
ScaleY = floor(Y);
close all;

ScaleCheck = 0;
disp('Use the cursor to select the scale range (lower & upper extremes).')

while ScaleCheck == 0
    % Set Distance Conversion Factor
    DistScale = OneFrame(ScaleY(1,1):ScaleY(2,1),ScaleX(1,1):ScaleX(2,1));
    DistCon = FrameCalibrate(DistScale);
    
    % Calibration Check
    CalCheck = input('Do you need to repeat the calibration (Y/N)? ','s');
    
    if lower(CalCheck) == 'n'
        ScaleCheck = 1;
    end
end

%% Calculate the Vessel Diameter
DiameterPixel = zeros(nFrames,1);

for k = 1:nFrames
    
    %Read Image at start of Video
    Vessel = read(USObj,k);
    
    if k == 1
        %Select ROI around the vessel 
        figure; imagesc(Vessel);
        [DROIx,DROIy,Center,theta] = VesselROI(Vessel);
    end
    
    % Restrict image area to ROI
    DiameterImage = Vessel(DROIy(1,1):DROIy(2,1),...
        DROIx(1,1):DROIx(2,1));
        
    %Calculate the Diameter
    DiameterPixel(k,1) = AutoDiameter(DiameterImage,Center',theta(1));
end

% Replace Gaps and Errant Values with Mean Diameter
MeanD = mean(DiameterPixel(intersect(find(DiameterPixel~= 200),find(DiameterPixel~=0))));
DiffD = diff(DiameterPixel(find(DiameterPixel~=0)))/MeanD;

for i = 1:length(DiameterPixel)-1
    if DiameterPixel(i,1) == 200
        DiameterPixel(i,1) = MeanD;
    end
end

% If Diameter Varies by > 20% from the Mean, Replace with Mean Value
for i = 1:length(DiffD)
    if DiffD(i,1) > 0.2
        DiameterPixel(i,1) = MeanD;
    end
end

% Convert Diameter to Centimeters
Diameter = DiameterPixel.*DistCon;
MeanDiam = mean(Diameter);
StdDiam = std(Diameter);

% Diameter Time
dt = 1/VidSamRate;
DTime = (0:dt:nFrames*dt-dt)';

% Filter Diameter - Smooth with Savitsky Golay Filter
DFilt = sgolayfilt(Diameter,3,11);
    
%% Plot and Summarize

% Error Limits - +/- 10% from Mean Diameter
% User Inspection of Data - If diameter exceeds error limits, repeat
% epoch selection or enter a single diameter value.
for i = 1:length(DFilt)
    PosError(i,1) = mean(DFilt) + .1*mean(DFilt);
    NegError(i,1) = mean(DFilt) - .1*mean(DFilt);
end

figure, grid on, hold on;
plot(DTime, DFilt,'r','DisplayName','Filtered Diameter', 'LineWidth', 2); hold on;
plot(DTime,PosError,'k','LineWidth',2); hold on;
plot(DTime,NegError,'k', 'LineWidth',2); hold off;
title('Vessel Diameter & +/- 10% Error Limits');
xlabel('Time(s)');
ylabel('Diameter (cm)');
legend('Diameter');
pause;

DistCon
MeanDiam
