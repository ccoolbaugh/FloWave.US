%% FloWaveUS - Ultrasound Blood Flow Analysis
%
% Purpose: Analysis of duplex ultrasound blood flow data from digital 
%          screen captures.
%
% Inputs: Digital Video of Ultrasound Data - AVI (PC) or MOV (MAC) formats
% Outputs: Blood flow time series and cardiac cycle data to csv files.
%
% Functions: FrameCalibrate, GapInterpolate, VesselROI, AutoDiameter, 
%            imgaussian (written by Dirk-Jan Kroon), ginputc (written by
%            Jiro Doke) - See associated license agreements for copyright
%            information, FigureLoopingFcn
%
% Help: See Documentation on Github for instruction manuals and example
% videos.
%
% version: 0.3.0
%
% Crystal Coolbaugh
% July 30, 2015
% Copyright 2015-2016 Crystal Coolbaugh
%% Format Workspace
clc
clear
close all;
format compact;
disp('PLEASE FOLLOW THE ON-SCREEN PROMPTS TO PROCESS THE ULTRASOUND DATA.')

%% Platform Calibration
% Users should create a platform calibration file for each experimental
% condition. The calibration ".mat" file can be created with
% PlatformCalibration.m. 

PlatSetFile = input('Enter the ultrasound setup filename (e.g. Settings.mat): ', 's');
load(PlatSetFile)

%% Identify Folder with Video
% Add video folder to the current path
disp('Choose the file directory that contains the US video.');
DirName = uigetdir; % Open browser to identify directory path
addpath(DirName);
ls(DirName) % List file contents for easy copy-paste of filename

%% Import Digital Video Data
% The digital video format depends on the video compression settings
% available in MATLAB. PC systems can use .avi or .mov formats. MAC systems
% will need to use .mov formats.

[USObj,nFrames,VidSamRate,Duration,OneFrame,BaseY1,FramesPerPeriod] = VideoData(PlatSetFile);

%% Pixel Calibration (FrameCalibrate Function)
% Determine the factor needed to convert pixels to the quantity of
% interest: velocity, time, and distance. 
% Position of the scales may need to be adjusted for different US
% equipment. Inspect a single image frame and adjust the zoomed portion.

[VelCon,TimeCon,DistCon] = PixelCal(OneFrame,PlatSetFile);

%% Identify Video Start Frame
% The ultrasound video should be edited with an external program so that
% the pulse wave data on the initial start frame of the video represents a
% complete update of data. To account for errors with editing, the user is
% presented an animation of the start of the video with the individual
% frame numbers. The user can press the  spacebar to stop the animation at
% the correct start frame. 
 
[StartFrame] = StartFrame(USObj,FramesPerPeriod);

%% Pulse Wave Data Extraction 
% For all frames in the video, this code a) defines an ROI for the pulse
% wave data, b) extracts the time average mean data (calculated teal line),
% c) extracts the max doppler envelope (TAMax), and e) scales 
% the velocity and time data from pixels to the appropriate units.

% NOTE: If time average mean data is not displayed in real-time, the mean 
% of the doppler envelope may be used as an alternative. 

[VelType,VelAll,MeanVelAll,WtMeanVelAll,EnvelopeVelAll] = PulseExtraction(StartFrame,FramesPerPeriod,nFrames,USObj,BaseY1,VelCon,PlatSetFile);

%% Create Composite Pulse-Wave Data
% Calls GapInterpolate function - replaces errant marker ("200") with
% linear interpolated velocity values.

[TimeFinal,VelInterp] = CompositePWData(TimeCon,VelAll,MeanVelAll,EnvelopeVelAll,WtMeanVelAll,VelType);

%% Define Data Epochs (VesselROI, AutoDiameter, imgaussian functions)
% The user interactively selects the end of data epochs (e.g. resting,
% contraction, post-contraction, recovery). Epochs should define shifts in
% the data frequency of waveform shape. Periods of movement in the
% ultrasound probe can also be identified for removal. 

[EpochTime,DRoiTime,DFilt,DTime] = DefineEpochs(TimeFinal,VelInterp,TimeCon,VidSamRate,nFrames,USObj,DistCon);

%% Create Common Time Series - Velocity & Diameter / Calculate Blood Flow and Shear Rate
% Diameter measures and time values are inherently sampled at different
% rates; therefore, a linear interpolation was applied to upsample the data
% to a common rate of 100 Hz. 

[FsNew,Time100,Vel100,Diam100,BloodFlow,Shear,BFFilt,SRFilt] = CommonData(VidSamRate,TimeCon,DTime,DFilt,TimeFinal,VelInterp);

%% Cardiac Cycle Analysis & Create Composite Data Set for All Epochs
% Identify the start/end of each cardiac cycle to determine the blood flow
% at peak systole, diastole, and end of diastole. Use the peak
% antegrade and retrograde shear to calculate the oscillatory shear index. 

[AllPSVTime,AllPSV,AllPDVTime,AllPDV,AllEDVTime,AllEDV,AllMBF1,AllData,AllTimeSeries,AnalysisSet] = AnalyzeEpochs(EpochTime,Time100,DRoiTime,BFFilt,SRFilt,Vel100,Diam100,FsNew,VelCon,TimeCon,DistCon,BaseY1,BloodFlow,Shear);

% Plot Cardiac Cycle Points and Mean Blood Flow
FinalPlots(AllPSVTime,AllPSV,AllPDVTime,AllPDV,AllEDVTime,AllEDV,AllMBF1);

%% Write to File
% Write processed cyclic data to file
% Write processed time series data to file
% Write analysis settings to file

WriteFile(AllData,AllTimeSeries,AnalysisSet);
