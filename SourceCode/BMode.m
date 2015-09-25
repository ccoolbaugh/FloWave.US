%% BMode
%
% Purpose: This code was written to analyze ultrasound (US) B-mode images
% from an analog video recording. The code uses an automated diameter 
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
% US Equipment: GE Logiq Book e
% Video Equipment: Sony Multi-Functional DVD Recorder (VRD-MC6)
% Video Conversion: VBO to AVI (PC) or MOV (MAC)
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

%% Platform Calibration
% If needed, the user can calibrate the program to the ultrasound screen
% dimensions. The calibration will identify ROIs for the velocity, time,
% and distance scales and the position of the pulse wave data. Selection of
% the color profile for the zero velocity position and the TAMean are also
% identified. 

PlatCal = menu('Use default platform calibration settings?', 'Yes', 'No');

if PlatCal == 1
    load GELogiqBookeDefault.mat
    disp('Loaded default settings')
else
    PlatSetFile = input('Enter the platform calibration settings filename (e.g. Settings.mat)','s');
    load PlatSetFile
    disp('Loaded custom settings')
end
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

OneFrame = read(USObj, 1);

ScaleCheck = 0;
disp('Use the cursor to select the scale range (lower & upper extremes).')

while ScaleCheck == 0
    % Set Distance Conversion Factor
    DistScale = OneFrame(DistY(1):DistY(2),DistX(1):DistX(2),:);
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
        [DROIx,DROIy,Center,Mask,theta] = VesselROI(Vessel);
    end
    
    % Restrict image area to ROI
    DiameterImage = Vessel(DROIy(1,1):DROIy(2,1),...
        DROIx(1,1):DROIx(2,1));
    
    % Mask Doppler gate
    DiameterImage((Center(1)-round((Mask/2))):(Center(1)+round((Mask/2))),:)=0;
    
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
Mask
theta
MeanDiam