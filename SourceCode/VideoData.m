function [USObj,nFrames,VidSamRate,Duration,OneFrame,BaseY1,FramesPerPeriod] = VideoData(PlatSetFile)
%%% Import Digital Video Data
% The digital video format depends on the video compression settings
% available in MATLAB. PC systems can use .avi or .mov formats. MAC systems
% will need to use .mov formats. 

load(PlatSetFile)

% Enter the Video Filename
VideoName = input('Type the video filename and extension (e.g. USvideo.avi): ','s');

% Construct a reader object using VideoReader function
USObj = VideoReader(VideoName);

% Set Image Frame Size - NOTE: Digital video resolution should match the
% platform calibration resolution to ensure proper scaling of pixel data.
VidHeight = USObj.Height;
VidWidth = USObj.Width;

% Total Number of Image Frames
nFrames = USObj.NumberofFrames;

% Video Sampling Rate
VidSamRate = USObj.FrameRate;

% Video End Time
Duration = USObj.Duration;

% Display a Single Video Frame for Visual Inspection of Doppler sweep speed,
% Velocity, Time, and Distance scales
OneFrame = read(USObj,1);
image(OneFrame)
axis image
title('Identify Doppler display sweep speed (max time range) & Calibration Scales', 'FontSize', 16);
disp('Inspect the image to identify the time scale, velocity scale, and diameter scale increments');

% Set the # of Frames per Screen Update for the US - Varies with Sweep Speed  
TimeMax = input('Enter the time (6 seconds = 6) needed to update the pulse wave data (e.g. the largest number on the time scale): ');
FramesPerPeriod = round(VidSamRate * TimeMax);

% Set Velocity Baseline - Zero Velocity Pixel Position
% [RGB] settings may need to be adjusted for different equipment. Use
% PlatformCalibration.m to determine the system's settings.
ROI(:,:,:) = OneFrame(PWy(1):PWy(2), PWx(1):PWx(2),:);

ROI_R = ROI(:,:,1);
ROI_G = ROI(:,:,2);
ROI_B = ROI(:,:,3);

CollapseR = mean(ROI_R,2);
CollapseG = mean(ROI_G,2);
CollapseB = mean(ROI_B,2);

MeanColorCollapse = [CollapseR CollapseG CollapseB];

BestMag = 100;
for k=1:length(CollapseR)
    vec_diff = BaseLineColor-MeanColorCollapse(k,:);
    mag_diff = sqrt(vec_diff(1)^2+vec_diff(2)^2+vec_diff(3)^2);
    if mag_diff<BestMag
        BaseY1 = k;
        BestMag=mag_diff;
    end;
end;


end

