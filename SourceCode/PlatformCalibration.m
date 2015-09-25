%% Platform Calibration - Ultrasound Calibration for FloWaveUS program
%
% Purpose: This code creates calibration settings for the postion of the
% ultrasound scales, the color of the pulse wave baseline, and the color of
% the pulse wave time average mean to use in the FloWaveUS program. 
%
% Inputs: Analog Video of Ultrasound Data - AVI (PC) or MOV (MAC) formats
% Outputs: MATLAB File containing calibration data for the ultrasound 
%          equipment 
%
% Crystal Coolbaugh
% July 28, 2015
% Copyright 2015 Crystal Coolbaugh

%% Format Workspace
clc
clear
close all;
format compact;
disp('PLEASE FOLLOW THE ON-SCREEN PROMPTS TO PROCESS THE ULTRASOUND DATA.')

%% Identify Folder with Video
% Add video folder to the current path
disp('Choose the file directory that contains the US video.');
DirName = uigetdir; % Open browser to identify directory path
addpath(DirName);
ls(DirName) % List file contents for easy copy-paste of filename

%% Import AVI Video Data
% Import an avi video file of the analog US data. 

% Enter the Video Filename
VideoName = input('Type the video filename and extension (e.g. USvideo.avi): ','s');

% Construct a reader object using VideoReader function
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

%% Set Scale Position ROIs

FrameOne = read(USObj,1);
CalScale = 0;

while ~CalScale
    figure;
    image(FrameOne);
    title('Define the Scale ROI: Velocity(Red), Time(Green), Distance(Yellow)')
    
    disp('Define the Velocity (Red), Time (Green), and Distance (Yellow) Scale ROI')
    disp('Click on 1. the upper left and 2. the lower right corners')
    
    %Red Cursor Select
    [VelX,VelY] =  ginputc(2,'Color','r', 'LineWidth',2);
    VelX = floor(VelX);
    VelY = floor(VelY);
    
    %Green Cursor Select
    [TimeX,TimeY] = ginputc(2,'Color','g', 'LineWidth',2);
    TimeX = floor(TimeX);
    TimeY = floor(TimeY);
    
    %Yellow Cursor Select
    [DistX,DistY] = ginputc(2,'Color','y', 'LineWidth',2);
    DistX = floor(DistX);
    DistY = floor(DistY);
    
    close all;

    %Review ROI Positions for Accuracy
    figure;
    image(FrameOne(VelY(1):VelY(2),VelX(1):VelX(2),:));
    title('Velocity Scale ROI');
    pause;
    close all;
    
    figure;
    image(FrameOne(TimeY(1):TimeY(2),TimeX(1):TimeX(2),:));
    title('Time Scale ROI');
    pause;
    close all;
    
    figure;
    image(FrameOne(DistY(1):DistY(2),DistX(1):DistX(2),:));
    title('Distance Scale ROI');
    pause;
    close all;
    
    %ROI Check
    ROITest = input('Were the scale positions selected correctly(Y/N)?','s');
    if lower(ROITest) == 'y'
        CalScale = 1;
    end
end

%% Set Pulse Wave ROI Size

PulseWave = 0;

while ~PulseWave
    figure;
    image(FrameOne);
    
    disp('Define the Pulse Wave Data (magenta) ROI');
    disp('Click on 1. the upper left and 2. the lower right corners of the ROI');
    disp('The ROI should include the full dynamic range of the Pulse Wave Data');
    
    %Magenta Cursor Select
    [PWx,PWy] = ginputc(2,'Color','m','LineWidth',2);
    PWx = floor(PWx);
    PWy = floor(PWy);
    PWHeight = diff(PWy);   %Define vector size for pulse wave data
    PWLength = diff(PWx);   
    
    close all;
    
    %Review ROI Position for Accuracy
    figure;
    PulseWaveData = FrameOne(PWy(1):PWy(2),PWx(1):PWx(2),:);
    image(PulseWaveData);
    title('Pulse Wave Data ROI');
    pause;
    close all;
    
    %ROI Check
    ROITest2 = input('Was the pulse wave data selected correctly (Y/N)?','s');
    if lower(ROITest2) == 'y'
        PulseWave = 1;
    end
end

%% Define Zero Axis for Pulse Wave Data

PulseZero = 0;

while ~PulseZero
    %Display Pulse Wave Data ROI
    figure;
    image(PulseWaveData);
    title('Pulse Wave Data ROI');
    
    %Select Position of Zero Velocity 
    disp('Click on 5 points along the Zero Velocity Position');
    [ZeroX,ZeroY] = ginputc(5,'Color','c','LineWidth',2);
    ZeroX = floor(ZeroX);
    ZeroY = floor(ZeroY);
    
    close all;
    
    %Average RGB values at each X position
    for i = 1:5
        Zero_R(i) = PulseWaveData(ZeroY(i),ZeroX(i),1);
        Zero_G(i) = PulseWaveData(ZeroY(i),ZeroX(i),2);
        Zero_B(i) = PulseWaveData(ZeroY(i),ZeroX(i),3);
    end
    
    BaseLine_R = mean(Zero_R);
    BaseLine_G = mean(Zero_G);
    BaseLine_B = mean(Zero_B);
    
    BaseLineColor = [BaseLine_R,BaseLine_G,BaseLine_B]
    
    %Baseline Check
    ColorTest = input('Was the baseline color selected correctly (Y/N)?','s');
    if lower(ColorTest) == 'y'
        PulseZero = 1;
    end
end

%% Define TAMean RGB Ratios

%Display TAMean extraction with defaul settings
TealLow = 0.25;
TealHigh = 0.7;

%Recast pulse wave data as type double 
ROI_db = double(PulseWaveData);
ROIsize = size(ROI_db);

%Extract TAMean pixel intensity
for i = 1:ROIsize(1)
    for j = 1:ROIsize(2)
        if ((ROI_db(i,j,1)/ROI_db(i,j,2) <= TealHigh) && ...
                (ROI_db(i,j,1)/ROI_db(i,j,2) >= TealLow) && ...
                (ROI_db(i,j,1)/ROI_db(i,j,3) <= TealHigh) && ...
                (ROI_db(i,j,1)/ROI_db(i,j,3) >= TealLow))
            
            I = 1/3*(ROI_db(i,j,1)+ROI_db(i,j,2)+ROI_db(i,j,3));
            TAMeanData(i,j,:) = I;
        end
    end
end

%Display TAMean Default 
figure;
image(TAMeanData);
title('Default TAMean Settings');

TAMeanReset = input('Do you want to reset the TAMean extraction ratios (Y/N)?', 's');

%Reset TAMean Ratios
if lower(TAMeanReset) == 'y'
    TACheck = 0;
    clear TAMeanData;
else
    TACheck = 1;
end

while ~TACheck
    
    %Display Pulse Wave Data ROI
    figure;
    image(PulseWaveData);
    title('Click on 10 points to define the TAMean Color');

    %Select Points to Define TA Mean Color
    disp('Click on 10 points to define the TAMean Color');
    [TAx,TAy] = ginputc(10,'Color','r','LineWidth',2);
    TAx = floor(TAx);
    TAy = floor(TAy);
    close all;
  
    %Identify RGB values at each position
    TA_R = zeros(1,10);
    TA_G = zeros(1,10);
    TA_B = zeros(1,10);
    RGRatio = zeros(1,10);
    RBRatio = zeros(1,10);
    
    for i  = 1:10
        TA_R(1,i) = ROI_db(TAy(i),TAx(i),1);
        TA_G(1,i) = ROI_db(TAy(i),TAx(i),2);
        TA_B(1,i) = ROI_db(TAy(i),TAx(i),3);
    end
    
    %Calculate R/G and R/B Ratios
    RGRatio = TA_R./TA_G;
    RBRatio = TA_R./TA_B;
    
    %Upper and Lower Ratio Range
    TealLow_New = min(min(RGRatio,RBRatio));
    TealHigh_New = max(max(RGRatio,RBRatio));
    
    %Display New TAMean Extraction
    for i = 1:ROIsize(1)
        for j = 1:ROIsize(2)
            if ((ROI_db(i,j,1)/ROI_db(i,j,2) <= TealHigh) && ...
                    (ROI_db(i,j,1)/ROI_db(i,j,2) >= TealLow) && ...
                    (ROI_db(i,j,1)/ROI_db(i,j,3) <= TealHigh) && ...
                    (ROI_db(i,j,1)/ROI_db(i,j,3) >= TealLow))
                
                I = 1/3*(ROI_db(i,j,1)+ROI_db(i,j,2)+ROI_db(i,j,3));
                TAMeanData(i,j,:) = I;
            end
        end
    end
    
    %Display TAMean Extraction
    figure;
    image(TAMeanData);
    title('New TAMean Settings');
    
    %Reselect or Use Defaults
    TAMeanReset2 = menu('Select TAMean extraction settings:', 'New', 'Default', 'Repeat selection');
    
    if TAMeanReset2 == 1
        TealLow = TealLow_New;
        TealHigh = TealHigh_New;
        TACheck = 1;
        close all;
    elseif TAMeanReset2 == 2
        TealLow = 0.25;
        TealHigh = 0.7
        TACheck = 1;
        close all;
    elseif TAMeanReset2 == 3
        disp('Repeat selection of TAMean color values.')
        clear TAMeanData;
        close all;
    end
end

%% Save Platform Calibration Settings

PlatSetFile = input('Enter a filename for the platform settings (e.g. Settings.mat): ', 's');

save(PlatSetFile, 'VelX','VelY','TimeX','TimeY', ...
    'DistX','DistY','PWx','PWy','PWHeight','PWLength','BaseLineColor',...
    'TealLow','TealHigh');


    
    
