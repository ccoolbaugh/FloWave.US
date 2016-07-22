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
% version: 0.2.1
%
% Crystal Coolbaugh
% July 30, 2015
% Copyright 2015 Crystal Coolbaugh
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

%% Pixel Calibration (FrameCalibrate Function)
% Determine the factor needed to convert pixels to the quantity of
% interest: velocity, time, and distance. 
% Position of the scales may need to be adjusted for different US
% equipment. Inspect a single image frame and adjust the zoomed portion.

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
%% Identify Video Start Frame
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

%% Pulse Wave Data Extraction 
% For all frames in the video, this code a) defines an ROI for the pulse
% wave data, b) extracts the time average mean data (calculated teal line),
% c) extracts the max doppler envelope (TAMax), and e) scales 
% the velocity and time data from pixels to the appropriate units.

% NOTE: If time average mean data is not displayed in real-time, the mean 
% of the doppler envelope may be used as an alternative. 

% Create Loop for Batch Processing
count = 1; % Count the # of Screen Updates

% Select Velocity Extraction Method - TAMean or Raw
VelType = menu('Extract the calculated TAMean or the raw Doppler spectrum?', 'TAMean', 'Raw');

for k = StartFrame:FramesPerPeriod:nFrames
    
    % Restrict Video Data to Pulse Wave ROI - defined in platform
    % calibration
    vidframe = read(USObj,k);
    ROI(:,:,:) = vidframe(PWy(1):PWy(2),PWx(1):PWx(2),:);
    
    % Allocate Storage for Velocity Data
    VelValue = zeros(1,PWLength);
    MaxI = zeros(1,PWLength);
    VelPos = zeros(1,PWLength);
    
    if VelType == 1
        %Recast ROI as type double
        ROI_db = double(ROI);
        
        % Extract Time Average Mean Data
        % Teal Pixels defined by Red/Green & Red/Blue Ratio:
        % Default Settings: Low = 0.25 & High = 0.7
        % [RGB] of the TAMean line may need to be adjusted for different
        % equipment - Use PlatformCalibration.m to determine the system
        % settings.
        
        % Allocate memory for data extraction
        TAMeanData = zeros(PWHeight,PWLength,1);
        
        % For all pixels in ROI, keep teal pixels and convert to intensity
        for i = 1:PWHeight
            for j = 1:PWLength
                if ((ROI_db(i,j,1)/ROI_db(i,j,2) <= TealHigh) && ...
                        (ROI_db(i,j,1)/ROI_db(i,j,2) >= TealLow) && ...
                        (ROI_db(i,j,1)/ROI_db(i,j,3) <= TealHigh) && ...
                        (ROI_db(i,j,1)/ROI_db(i,j,3) >= TealLow))
                    
                    I = 1/3*(ROI_db(i,j,1)+ROI_db(i,j,2)+ROI_db(i,j,3));
                    TAMeanData(i,j,:) = I;
                end
            end
        end
        
        % Find maximum intensity within TAMean Data and convert position to a
        % [time, velocity] data point
        for i = 1:PWLength
            % Find Max Intensity per Column
            MaxI(1,i) = max(TAMeanData(:,i));
            
            if MaxI(1,i) > 30   %Default threshold for Intensity = 30; Lower intensity values add noise
                VelIndex = find(TAMeanData(:,i) == MaxI(1,i));
                
                VelPos(1,i) = median(VelIndex); % Median accounts for multiple maximum
                VelValue(1,i) = (BaseY1 - VelPos(1,i))*VelCon;
            else
                VelValue(1,i) = 200;    % Mark missing data for interpolation
            end
        end

         % Add Data to Array
         VelAll(:,count) = VelValue';
         
    else
        % Extract Doppler Envelope Data - grayscale threshold (GrayThresh
        % set in PlatformCalibration.m)
        % For all pixels in ROI, keep gray pixels
        RawData = zeros(PWHeight,PWLength,1);
                
        % Mask velocity baseline to prevent noise in velocity envelope
        BaseMask = 2;
        ROI((BaseY1-BaseMask:BaseY1+BaseMask),:,:) = 0;
        
        %Recast ROI as type double
        ROI_db = double(ROI);
        
        % Threshold for Gray Pixels (RGB must be within threshold)
        for i = 1:PWHeight
            for j = 1:PWLength
                if (abs(ROI_db(i,j,1) - ROI_db(i,j,2)) <= GrayThresh) && (abs(ROI_db(i,j,2)...
                        -ROI_db(i,j,3)) <= GrayThresh)
                    I = 1/3*(ROI_db(i,j,1)+ROI_db(i,j,2)+ROI_db(i,j,3));
                    
                    RawData(i,j,:) = I;
                end
            end
        end
        
        % Assign Velocity Pixel Positions
        for i = 1:PWLength
            % Filter Noisy Pixels
            IntIndex = find(RawData(:,i) >= 10);
            
            %Calculated intensity weighted mean and mean velocity values.
            %Weighted mean is currently an exploratory measure. 
            if IntIndex ~= 0
                % Calculate Intensity Weighted Mean Position
                PixelVec = (1:length(IntIndex))';
                VelVec = RawData(IntIndex,i);
                WtMeanPixel = round(sum(PixelVec.*VelVec)./sum(VelVec));
                WtMeanPos(1,i) = IntIndex(WtMeanPixel);
                
                % Calculate Intensity Weighted Mean Velocity
                WtMeanVelValue(1,i) = (BaseY1 - WtMeanPos(1,i))*VelCon;
                
                % Calculate Velocity Envelope
                for j = 1:length(IntIndex)
                    VelIndex(j,1) = (BaseY1 - IntIndex(j,1))*VelCon;
                end
                
                % Assign Max/Min/Mean Velocities
                MaxVelValue(1,i) = max(VelIndex);
                MinVelValue(1,i) = min(VelIndex);
                MeanVelValue(1,i) = mean(VelIndex);
            else
                % Mark missing data for interpolation
                WtMeanVelValue(1,i) = 0;
                MeanVelValue(1,i) = 0;
                MaxVelValue(1,i) = 200;
                MinVelValue(1,i) = -200;
            end
            
            %Calculate Pulse Wave Spectrum Envelope
            EnvelopeVel(1,i) = MaxVelValue(1,i) + MinVelValue(1,i);
            
            % Reset Index 
            clear IntIndex VelIndex
        end
        
         % Add Data to Array
         MeanVelAll(:,count) = MeanVelValue';
         WtMeanVelAll(:,count) = WtMeanVelValue';
         EnvelopeVelAll(:,count) = EnvelopeVel';
    end
    
    % Increment Total Frame Count
    TotalFrameCount = count;
    count = count + 1;
    
    % Clear variables to avoid overlap with next image
    clear TAMeanData RawData I 
 
end

%% Create Composite Pulse-Wave Data
% Calls GapInterpolate function - replaces errant marker ("200") with
% linear interpolated velocity values.

% Reshape Velocity Vectors & Correct Missing Data Gaps
if VelType == 1
    TimeFinal = (0:TimeCon:numel(VelAll)*TimeCon-TimeCon);  %Time Vector
    VelFinal = reshape(VelAll,numel(VelAll),1);
    VelInterp = GapInterpolate(VelFinal,TimeFinal,TimeCon);
else
    TimeFinal = (0:TimeCon:numel(MeanVelAll)*TimeCon-TimeCon); % Time Vector
    MeanVelFinal = reshape(MeanVelAll,numel(MeanVelAll),1);
    VelInterp = GapInterpolate(MeanVelFinal,TimeFinal,TimeCon);
    EnvelopeVelFinal = reshape(EnvelopeVelAll,numel(EnvelopeVelAll),1);
    MaxVelInterp = GapInterpolate(EnvelopeVelFinal,TimeFinal,TimeCon);
    
    % Calculate Intensity-Weighted Mean (Exploratory Analysis)
    WtMeanVelFinal = reshape(WtMeanVelAll,numel(WtMeanVelAll),1);
    WtMeanVelInterp = GapInterpolate(WtMeanVelFinal,TimeFinal,TimeCon);
end

%% Define Data Epochs (VesselROI, AutoDiameter, imgaussian functions)
% The user interactively selects the end of data epochs (e.g. resting,
% contraction, post-contraction, recovery). Epochs should define shifts in
% the data frequency of waveform shape. Periods of movement in the
% ultrasound probe can also be identified for removal. 

% Identify the Epoch Boundaries
disp('Click on the pulse wave data to define sections (epochs) of the experiment. Press enter when finished.')
EpochCheck = 0;
ROICheck = 0;

while ROICheck == 0
    while EpochCheck == 0
        % Select Data Epochs
        figure, plot(TimeFinal, VelInterp, 'k');
        xlabel('Time (s)'); ylabel('Velocity (cm/s)');
        title('Use the Cursor to Select Epochs. Press enter when done.')
        
        [t_epoch, v_epoch] = ginput;
        T = TimeCon;
        close all;
        
        disp('Number of epochs selected: ')
        disp(length(v_epoch))
        
        for k = 1:length(t_epoch)
            if k == 1
                ELow = 1;
                EHigh = intersect(find(TimeFinal<=(t_epoch(k)+T)),find(TimeFinal>=t_epoch(k)-T));
            else
                % Find time index for lesser and greater time points
                ELow = intersect(find(TimeFinal<=(t_epoch(k-1)+T)),find(TimeFinal>=t_epoch(k-1)-T));
                EHigh = intersect(find(TimeFinal<=t_epoch(k)+T),find(TimeFinal>=t_epoch(k)-T));
            end
            
            if ~isempty(EHigh)
                EMean = VelInterp(ELow(1):EHigh(1));
                ETime = TimeFinal(ELow(1):EHigh(1));
                
                % Display Selected Epoch
                figure;
                plot(ETime, EMean, 'k');
                xlabel('Time (sec)');
                ylabel('Velocity (cm/s)');
                title(['Data Selected for Epoch: ' int2str(k)]);
                pause;
                
                % Assign Data to Cell Arrays
                EpochTime{k} = ETime;
                
                % Identify Time Points for Diameter ROI Selection
                % Correct for delay in Velocity/Number of Frames and
                % diameter display
                DRoiTime(k,1) = ETime(end);
                
            else
                disp('ERROR: Final point selection is outside of the data limits. Reselect the Epochs.');
            end
        end
        
        % Visual Error Inspection of Epoch Selection
        EpochTest = input('Were Epochs Selected Correctly (Y/N)?', 's');
        if lower(EpochTest) == 'y'
            disp('Advancing to Diameter Selection')
            close all;
            EpochCheck = 1;
        else
            disp('Reselect Epochs.')
            clear t_epoch v_epoch
        end
    end
    
    % VESSEL DIAMETER MEASUREMENT
    % Options: Use a single value or use automated algorithm
    % Single Value: User enters a single diameter value 
    % Automated: The user is shown an image of the vessel at the start of each data
    % epoch. A rectangular ROI is selected to include the near and far
    % vessel walls. The user then defines the center of the vessel and the
    % vessel angle. 
    
    DiameterOption = menu('Choose a vessel diameter method:', 'Single Value', 'Automated');
    
    if DiameterOption == 1
        % Create Diameter Time Vector
        dt = 1/VidSamRate;
        DiamEndFrame = nFrames; 
        DTime = (0:dt:DiamEndFrame*dt-dt)';
        
        % Enter diameter value
        DFilt = zeros(length(DTime),1);
        DFilt(:,1) = input('Enter a diameter value (cm) to use for blood flow calculations: ');
        
        ROICheck = 1;
        disp('Advancing to Cycle Analysis')
    elseif DiameterOption == 2
        
        % WARNING: Automated diameter is still a beta-version analysis
        disp('WARNING: Automated measurement of the vessel diameter is a BETA analysis.')
        
        % Create Fiducial Markers for Epochs - Define Vessel ROI at start of Each Epoch
        ROIMark = [1;round(DRoiTime*VidSamRate)];
        
        % Set Vector Length to Match the Velocity Vector
        DiamEndFrame = nFrames; 
        DiameterPixel = zeros(DiamEndFrame,1);
        
        for k = 1:DiamEndFrame
            % Read Image at start of Epoch
            Vessel = read(USObj,k);
            
            for j = 1:length(ROIMark)
                if k == ROIMark(j)
                    figure; imagesc(Vessel);
                    title(['Image Time: ' int2str(ROIMark(j)/VidSamRate)]);
                    ROIChoice = menu('Select an ROI?','Yes', 'No');
                    if ROIChoice == 1
                        Contraction = 0;
                        [DROIx,DROIy,Center,theta] = VesselROI(Vessel);
                    else
                        Contraction = 1; % No vessel measurement during contraction
                    end
                end
            end
                 
            if Contraction == 0
                % Restrict image to ROI
                DiameterImage = Vessel(DROIy(1,1):DROIy(2,1),...
                    DROIx(1,1):DROIx(2,1));
                % Calculate the Diamter
                DiameterPixel(k,1) = AutoDiameter(DiameterImage,Center',theta(1));
            else
                DiameterPixel(k,1) = 0;
            end
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
        
        % Create Diameter Time Vector
        dt = 1/VidSamRate;
        DTime = (0:dt:DiamEndFrame*dt-dt)';
        
        % Filter Diameter - Smooth with Savitsky Golay Filter
        DFilt = sgolayfilt(Diameter,3,11);
        
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
        
        % Visual Error Inspection of Diameter Detection
        ROITest = menu('Repeat Epoch Selection to Improve Diameter Detection?', 'Yes', 'No','Use Single Diameter Value');
        if ROITest == 1
            disp('Reselect Epochs.')
            EpochCheck = 0;
            clear t_epoch v_epoch ROIMark DRoiTime EpochTime
        elseif ROITest == 2
            disp('Advancing to Cycle Analysis')
            close all;
            ROICheck = 1;
        elseif ROITest == 3
            DFilt = zeros(length(DTime),1);
            DFilt(:,1) = input('Enter a diameter value (cm) to use for blood flow calculations: ');
            close all;
            ROICheck = 1;
            disp('Advancing to Cycle Analysis')
        end
    end
end

%% Create Common Time Series - Velocity & Diameter
% Diameter measures and time values are inherently sampled at different
% rates; therefore, a linear interpolation was applied to upsample the data
% to a common rate of 100 Hz. 

% Linear interpolate diameter and velocity time series to 100 Hz
FsNew = 100;    % New Sampling Rate
FsD = VidSamRate;   % Current Diameter Sampling Rate
FsV = 1/TimeCon;    % Current Velocity Sampling Rate

% Set Endpoint for New Time Vector - Removes possible artifact at end of
% diameter time series due to epoch selection
if DTime(end) < TimeFinal(end)
    TEnd = DTime(end);
else
    TEnd = TimeFinal(end);
end

% Create New Time Vector
Time100 = linspace(0,TEnd,TEnd*FsNew);

% Resample Diameter and Velocity Vectors
Vel100 = interp1(TimeFinal,VelInterp,Time100);
Diam100 = interp1(DTime,DFilt,Time100);
%% Calculate Blood Flow and Shear Rate

% Create Blood Flow (ml/min) and Shear Rate (1/s) Time Series
BloodFlow = zeros(length(Diam100),1);
Shear = zeros(length(Diam100),1);

for i = 1:length(Diam100)
    BloodFlow(i,1) = Vel100(i)*pi*(Diam100(i)/2)^2*60;  %ml/min
    if Diam100(i) ~= 0 
        Shear(i,1) = (4*Vel100(i))/Diam100(i);          %s^-1
    else
        Shear(i,1) = 0;
    end
end

% Smooth Time Series - Savitsky Golay Filter - remove noise 
BFFilt = sgolayfilt(BloodFlow,3,11);
SRFilt = sgolayfilt(Shear,3,11);

%% Cardiac Cycle Analysis
% Identify the start/end of each cardiac cycle to determine the blood flow
% at peak systole, diastole, and end of diastole. Use the peak
% antegrade and retrograde shear to calculate the oscillatory shear index. 

PkCountTotal = 0;
AnalyzeSuccess = 1; % Enable repeat analysis for each epoch

while AnalyzeSuccess <= length(EpochTime)
    
    if AnalyzeSuccess <= length(EpochTime)
        k = AnalyzeSuccess; % Increase counter if analysis completed correctly;
    else
        break
    end
    
    % Partition Filtered Blood Flow and Shear Data according to Epochs
    IndexE = find(Time100<=DRoiTime(k));
    EpochEnd(k,1) = IndexE(end);
    
    if k == 1
        EpochStart = 1;
    else
        EpochStart = EpochEnd((k-1),1);
    end
    
    BF = BFFilt(EpochStart:EpochEnd(k,1));
    SR = SRFilt(EpochStart:EpochEnd(k,1));
    Vel = Vel100(EpochStart:EpochEnd(k,1));
    Diam = Diam100(EpochStart:EpochEnd(k,1));
    OneTime = Time100(EpochStart:EpochEnd(k,1));
    
    % Plot Data Epoch - TAMean Velocity Data 
    figure;
    plot(OneTime, Vel, 'k');
    title(['Velocity vs. Time for Epoch: ' int2str(k)])
    
    % Choose to analyze or skip (muscle contraction or probe movement) data
    analyze = menu('Do you want to analyze this data?', 'Yes', 'No');
    close all;
    
    if analyze == 1
        
        % Smooth Artifacts for Peak Detect
        w = 16;     % Level of Smoothing 
        smooth = (ones(1,w) / w)';
        VelFilt1 = conv(Vel,smooth,'same');
                
        % Identify Systolic Peaks - Use Median Filtered Velocity Data
        PkCheck = 0;
        pkht = 20;  % Minimum vertical height for each peak; Default = 20 cm/s
        
        while PkCheck == 0;
            [peaks, locs] = findpeaks(VelFilt1,'minpeakheight', pkht);
            
            % Assess peak detection quality
            figure; 
            plot(VelFilt1,'k'); hold on;
            plot(locs,VelFilt1(locs), 'o','MarkerFaceColor','r', ... 
                'MarkerEdgeColor','r', 'MarkerSize', 5);
            title('Systolic Peak Detection')
            xlabel('Datapoint')
            ylabel('Velocity (cm/s)')
            hold off;
            
            QTest = input('Were the peaks correctly identified (Y/N)? ','s');
            
            if lower(QTest) == 'y'
                pkht_keep = pkht;
                pkht = 20;
                close all;
                PkCheck = 1;
            else
                pkht = input('Choose a new value for the minimum peak height (default = 20 cm/s): ');
            end
        end
        
        % Remove Double-Peaks (Artifacts)
        PkCheck2 = 0;
        pkwdth = 20; %Minimum distance between peaks; Default = 20 data points
        DblTest = locs; 
        
        while PkCheck2 == 0
            for i = 2:numel(DblTest)
                PkDiff = DblTest(i) - DblTest(i-1);
                if PkDiff <= pkwdth
                    DblTest(i) = NaN;
                    DblTest(i-1) = NaN;
                end
            end
            
            index = 1;
            for i = 1:numel(DblTest)
                if isnan(DblTest(i)) == 0
                    locs_new(index) = DblTest(i);
                    peaks_new(index) = peaks(i);
                    index = index+1;
                end
            end
            
            DblLocs = locs_new;
            DblPeaks = peaks_new;
            
            % Assess quality of peak distance threshold
            figure; 
            plot(VelFilt1,'k'); hold on;
            plot(DblLocs, VelFilt1(DblLocs),'o','MarkerFaceColor', 'r', ...
                'MarkerEdgeColor', 'r', 'MarkerSize',5);
            title('Double Peak Removal')
            xlabel('DataPoint');
            ylabel('Velocity (cm/s)')
            hold off;
            
            QTest2 = input('Were double-peaks correctly removed (Y/N)? ','s');
            
            if lower(QTest2) == 'y'
                pkwdth_keep = pkwdth;
                pkwdth = 20;
                close all;
                PkCheck2 = 1;
            else
                clear PeakDiff locs_new peaks_new   % Reset peak variables to avoid overlap
                DblTest = locs;
                pkwdth = input('Choose a new value for the minimum distance between peaks (default = 20 data points): ');
            end
        end
        
        % Find Peak Locations Filtered Data - Blood Flow and Shear Rate
        BFLocs = DblLocs;
        BFPeaks = BFFilt(DblLocs);
        
        SRLocs = DblLocs;
        SRPeaks = SRFilt(DblLocs);
        
        % DEFINE CARDIAC CYCLE FREQUENCY
        % Calculate time between peaks. 
        
        RRInt = diff(Time100(DblLocs));
        MeanRR = mean(RRInt);
        count = 1;
        
        % Exclude if +/- 30% different from mean RR
        for i = 1:length(RRInt)
            if (RRInt(i) <= (MeanRR+.3*MeanRR)) && (RRInt(i) >= (MeanRR-.3*MeanRR))
                RRIntKeep(count) = RRInt(i);
                count = count + 1;
            end
        end
        
        RRDur = mean(RRIntKeep);   % New mean for windowing 
            
        disp('Mean Cardiac Cycle Duration (seconds) is: ')
        disp(RRDur)
            
        % EXTRACT POINTS OF INTEREST FOR EACH CARDIAC CYCLE
        CycleCount = 1;
        
        % Allocate Memory
        Size = numel(DblPeaks); % # of Cycles to Analyze
        
        ISV = zeros(1,Size);    % Start of Systole
        ISVTime = zeros(1,Size);
       
        PSV = zeros(1,Size);    % Peak Systole
        PSVTime = zeros(1,Size);
        
        PDV = zeros(1,Size);    % Nadir of Diastole
        PDVTime = zeros(1,Size);
        
        EDV = zeros(1,Size);    % End of Diastole
        EDVTime = zeros(1,Size);
    
        Slope = zeros(1,Size);  % Systolic Slope
        SlopeTime = zeros(1,Size);
        
        MBF1 = zeros(1,Size); % Mean Blood Flow - Method 1
        MBF2 = zeros(1,Size);  % Mean Blood Flow - Method 2
        SysTime = zeros(1,Size);    % Time in Systole
        DiasTime = zeros(1,Size);   % Time in Diastole
        OSI = zeros(1,Size);    % Oscillatory Shear Index
     
        for i = 1:numel(DblPeaks)
            
            % Create Time Window
            WinDur = (RRDur)/2;
            
            Center = OneTime(DblLocs(i));
            Left = find(OneTime >= ((Center-WinDur) - 1/FsNew));
            Right = find(OneTime <= ((Center+WinDur) + 1/FsNew));
            
            WinTime = OneTime(Left(1):Right(end));
            WinBF = BFFilt(Left(1):Right(end));
            WinSR = SRFilt(Left(1):Right(end));
            WinVel = Vel(Left(1):Right(end));
            WinDiam = Diam(Left(1):Right(end));
            
            % Locate Peak Systole
            PSVLoc = find(WinTime == OneTime(DblLocs(i)));
            
            % Index timepoint with the waveform in the complete time series
            PSVt = find(Time100 == WinTime(PSVLoc));
            PSVTime(1,i) = Time100(PSVt);
            SlopeTime(1,i) = Time100(PSVt);
            
            % Locate Nadir Diastole - Limit to data following peak systole
            WinBFShort = -WinBF(PSVLoc:end);
            trht = 5;   % minimum vertical height for each peak (5 ml/min)
            [Trough,TLocShort] = findpeaks(WinBFShort,'minpeakheight', trht);
            
            % If multiple troughs are detected, keep the one that occurs
            % first in the time series
            if numel(TLocShort) > 1
                TLocShort = TLocShort(1);
            end
            
            % If no troughs are detected, assign a value of zero.
            if isempty(TLocShort) 
                PDVLoc = 0;
                
                % Index timepoint to match PSV
                PDVTime(1,i) = Time100(PSVt);
            else
                % Identify location in complete cardiac cycle
                PDVLoc = find(-WinBF == WinBFShort(TLocShort));
                
                % Index timepoint with the waveform in the complete time
                % series
                PDVt = find(Time100 == WinTime(PDVLoc));
                PDVTime(1,i) = Time100(PDVt);
            end
            
            % Locate start of systole
            PtDiff = diff(WinBF);
            
            % Assume start occurs within 0.25 s of Peak Systole 
            % (25 datapoints) or start of window
            if (PSVLoc - 25) >= 1
                for j = (PSVLoc-25):PSVLoc
                    if (PtDiff(j) >= 1.5) && ((j-1) >= 5)
                        Start = WinTime(j-1);
                        break
                    else
                        Start = WinTime(PSVLoc-20);
                    end
                end
            else
                for j = 1:PSVLoc
                    if (PtDiff(j) >= 1.5) && ((j-1) >= 5)
                        Start = WinTime(j-1);
                        break
                    else
                        Start = WinTime(1);
                    end
                end
            end
            
            ISVLoc = find(WinTime == Start);
            
            %If no start point is found, set to beginning of window
            if isempty(ISVLoc) == 1
                ISVLoc = 1;
            end
            
            % Index Timepoint with Waveform in Complete time Series
            ISVt = find(Time100 == WinTime(ISVLoc));
            ISVTime(1,i) = Time100(ISVt);
            
            % Locate end of window - placeholder for end of diastole
            EDVLoc = numel(WinTime);
            
            % Index timepoint with waveform in complete time series
            EDVt = find(Time100 == WinTime(EDVLoc));
            EDVTime(1,i) = Time100(EDVt);
            
            % Advance to Next Epoch
            AnalyzeSuccess = AnalyzeSuccess+1;
            
            %{
            % Quality Check - Three Cycles
            % If points are poorly selected, a new frequency should be
            % used and analysis repeated. 
            if (numel(DblPeaks) >= 3)
                if i <= 3
                    if PDVLoc ~= 0
                        KeyPts = [ISVLoc, PSVLoc, PDVLoc,EDVLoc];
                    else
                        KeyPts = [ISVLoc,PSVLoc,EDVLoc];
                    end
                    
                    figure;
                    plot(WinVel,'k'); hold on;
                    plot(KeyPts,WinVel(KeyPts),'o','MarkerFaceColor','r',...
                        'MarkerEdgeColor','r','MarkerSize',5);
                    title(['Data Extraction for Cycle: ' int2str(i)])
                    xlabel('Datapoint')
                    ylabel('Velocity (cm/s)')
                    hold off;
                    pause;
                    close all;
                elseif i == numel(DblPeaks)
                    Check = input('Were points correctly identified (Y/N)? ','s');
                    if lower(Check) == 'y'
                        AnalyzeSuccess = AnalyzeSuccess + 1;
                    else
                        break
                        disp('Repeat analysis with new peak find settings.')
                    end
                end
            else
                AnalyzeSuccess = AnalyzeSuccess+1;
            end
            %}
            
            % Store Points of Interest for Each Cycle
            ISV(1,i) = WinBF(ISVLoc);
            PSV(1,i) = WinBF(PSVLoc);
            if PDVLoc ~= 0
                PDV(1,i) = WinBF(PDVLoc);
            else
                PDV(1,i) = 0;
            end
            EDV(1,i) = WinBF(EDVLoc);
   
            % Calculate Slope and Oscillatory Shear Index
            Slope(1,i) = (WinBF(PSVLoc) - WinBF(ISVLoc))/(WinTime(PSVLoc) - WinTime(ISVLoc));
            OSI(1,i) = abs(min(WinSR))/(abs(max(WinSR))+abs(min(WinSR)));
        end
        
        % Assign End Diastole Values/Time Points
        % Based on Timepoints of the Initial Systole.
        ISVDiff = diff(ISVTime);    % Time between cardiac cycles
        RRLim = .3*RRDur;           % Cardiac Cycle Limit (30% of Mean RR)
        RRDurHigh = RRDur+RRLim; % Set High Cardiac Cyclic Limit
        RRDurLow = RRDur-RRLim;  %Set Low Cardiac Cycle Limit
        
        for i = 1:(numel(DblPeaks)-1)
            % Assign end of diastole as the next start of systole if
            % gap between cardiac cycle limits else assign NaN
            if (ISVDiff(i)<=RRDurHigh) && (ISVDiff(i)>= RRDurLow)
                EDVTime(1,i) = ISVTime(1,i+1);
                EDV(1,i) = ISV(1,i+1);
            else
                EDV(1,i) = NaN;
                EDVTime(1,i) = NaN;
            end
        end
        
        % Remove Gaps (NaNs) in the Data Set
        ISV = ISV(~isnan(EDV));
        ISVTime = ISVTime(~isnan(EDV));
        PSV = PSV(~isnan(EDV));
        PSVTime = PSVTime(~isnan(EDV));
        PDV = PDV(~isnan(EDV));
        PDVTime = PDVTime(~isnan(EDV));
        EDV = EDV(~isnan(EDV));
        EDVTime = EDVTime(~isnan(EDVTime));
        Slope = Slope(~isnan(EDV));
        SlopeTime = SlopeTime(~isnan(EDV));
        OSI = OSI(~isnan(EDV));
        
        % Calculate Mean BF weighted to time spent in systole & diastole
        for m = 1:numel(EDV)
            % Method 1: Mean of All BF from start to end of cardiac cycle
            % Method 2: Mean of BF weighted to time spent in systole & diastole
            % Systole (Exploratory Analysis)
            SysTime(1,m) = PSVTime(1,m) - ISVTime(1,m);
            SysStart = find(Time100 == ISVTime(1,m));
            SysEnd = find(Time100 == PSVTime(1,m));
            
            % Diastole
            DiasTime(1,m) = EDVTime(1,m) - PSVTime(1,m);
            DiasEnd = find(Time100 == EDVTime(1,m));
            
            % Cardiac Cycle Period
            Period = Time100(SysStart:DiasEnd);
            
            % Method 1:
            MBF1(1,m) = trapz(Period,BFFilt(SysStart:DiasEnd))*(1/(Time100(DiasEnd)-Time100(SysStart)));
            
            % Method 2: 
            MBF2(1,m) = SysTime(1,m)*(trapz(Period,BFFilt(SysStart:DiasEnd))*(1/(Time100(DiasEnd)-Time100(SysStart))))+ ...
                DiasTime(1,m)*(trapz(Period,BFFilt(SysStart:DiasEnd)*(1/(Time100(DiasEnd)-Time100(SysStart)))));
        end
        
        % Remove trailing zeros in SysTime, Diastime, MBF1, MBF2
        if (numel(EDV) < numel(DblPeaks))
            SysTime = SysTime(SysTime~=0);
            DiasTime = DiasTime(DiasTime~=0);
            MBF1 = MBF1(MBF1~=0);
            MBF2 = MBF2(MBF2~=0);
        end
        
        % Store Peak Find Settings
        PkHtSet(k) = pkht_keep;
        PkWdSet(k) = pkwdth_keep;
        
        % Store Cardiac Cycle Duration
        RRDurKeep(k) = RRDur;
        
        % Peak Amplitude and Count Per Epoch
        PeakCount(k) = length(DblPeaks);
        PkCountTotal = PkCountTotal + PeakCount(k);
        
        % Calibration Settings - Write only for first epoch
        if k == 1
            VelCal(k) = VelCon;
            TimeCal(k) = TimeCon;
            DistCal(k) = DistCon;
            ZeroCal(k) = BaseY1;
        else
            VelCal(k) = 0;
            TimeCal(k) = 0;
            DistCal(k) = 0;
            ZeroCal(k) = 0;
        end
        
    else
        
        % Assign Time for Epoch
        OneTime = Time100(EpochStart:(EpochEnd(k,1)-1));
        
        % Assign Zeroes for Each Variable
        for i = 1:numel(OneTime);
            PSVt = find(Time100 == OneTime(i));
            
            if numel(PSVt) > 1
                PSVt = PSVt(1);
            end
            
            ISV(1,i) = 0;
            ISVTime(1,i) = Time100(PSVt);
            PSV(1,i) = 0;
            PSVTime(1,i) = Time100(PSVt);
            PDV(1,i) = 0;
            PDVTime(1,i) = Time100(PSVt);
            EDV(1,i) = 0;
            EDVTime(1,i) = Time100(PSVt);
            Slope(1,i) = 0;
            SlopeTime(1,i) = Time100(PSVt);
            MBF1(1,i) = 0;
            MBF2(1,i) = 0;
            SysTime(1,i) = 0;
            DiasTime(1,i) = 0;
            OSI(1,i) = 0;
            RRIntKeep(1,i) = 0;
        end
        
         % Store Peak Find Settings
        PkHtSet(k) = 0;
        PkWdSet(k) = 0;
        
        % Store Cardiac Cycle Duration
        RRDurKeep(k) = 0;
        
        % Peak Amplitude and Count Per Epoch
        PeakCount(k) = 0;
        
        % Calibration Settings - Write only for first epoch
        if k == 1
            VelCal(k) = VelCon;
            TimeCal(k) = TimeCon;
            DistCal(k) = DistCon;
            ZeroCal(k) = BaseY1;
        else
            VelCal(k) = 0;
            TimeCal(k) = 0;
            DistCal(k) = 0;
            ZeroCal(k) = 0;
        end
        
        AnalyzeSuccess = AnalyzeSuccess+1;
    end
    
    % Store Data in Structures
    S(k).epoch = k;
    S(k).ISV = ISV;
    S(k).PSV = PSV;
    S(k).PDV = PDV;
    S(k).EDV = EDV;
    S(k).MBF1 = MBF1;
    S(k).MBF2 = MBF2;
    S(k).SysTime = SysTime;
    S(k).DiasTime = DiasTime;
    S(k).Slope = Slope;
    S(k).OSI = OSI;
    S(k).RRIntKeep = RRIntKeep;
    
    % Store time indices for points of interest
    S(k).ISVt = ISVTime;
    S(k).PSVt = PSVTime;
    S(k).PDVt = PDVTime;
    S(k).EDVt = EDVTime;
    S(k).Slopet = SlopeTime;
    S(k).Window = EDVTime - ISVTime;    
    
    % Clear Variables to Avoid Overlap with Next Epoch
    clear BF SR OneTime peaks locs ...
        RRIntKeep DblPeaks DblLocs Fd locs_new peaks_new vel

end

%% Create Composite Data Set for All Epochs

% Concatenate Data with the Structure
AllISV = [S(1,:).ISV];
AllISVTime = [S(1,:).ISVt];
AllPSV = [S(1,:).PSV];
AllPSVTime = [S(1,:).PSVt];
AllPDV = [S(1,:).PDV];
AllPDVTime = [S(1,:).PDVt];
AllEDV = [S(1,:).EDV];
AllEDVTime = [S(1,:).EDVt];
AllMBF1 = [S(1,:).MBF1];
AllMBF2 = [S(1,:).MBF2];
AllSysTime = [S(1,:).SysTime];
AllDiasTime = [S(1,:).DiasTime];
AllSlope = [S(1,:).Slope];
AllSlopeTime = [S(1,:).Slopet];
AllOSI = [S(1,:).OSI];
AllRRIntKeep = [S(1,:).RRIntKeep];
AllWindow = [S(1,:).Window];    

% Plot Cardiac Cycle Points 
figure, grid on, hold on;
plot(AllPSVTime,AllPSV, '--go','DisplayName', 'Peak Systole', 'LineWidth', 1.1, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 4); hold all;
plot(AllPDVTime,AllPDV,'--bo','DisplayName','Nadir Diastole','LineWidth',1.1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize', 4);hold all;
plot(AllEDVTime,AllEDV,'--mo','DisplayName','End Diastole','LineWidth',1.1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize', 4);hold off;
title('Blood Flow vs. Time','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Blood Flow (ml/min)','FontSize',14)
legend('Peak Systole','Diastole','End Diastole')

% Plot Mean Blood Flow
figure, axis tight; 
plot(AllPSVTime,AllMBF1,'k','DisplayName','Mean Blood Flow','LineWidth', 1.1, 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 4), grid on;
xlabel('Time (s)', 'FontSize', 14),
ylabel('Mean Blood Flow (ml/min', 'FontSize', 14),

%% Write to File

% Aggregate Cardiac Cycle Data  - Format for Validation Study Only
AllData = [AllPSVTime',AllPSV',AllPDVTime',AllPDV',AllEDVTime',AllEDV', ...
    AllISVTime',AllISV',AllMBF1',AllMBF2',...
    AllSysTime',AllDiasTime',AllOSI',AllWindow'];
 
% Aggregate Time Series Data
AllTimeSeries = [Time100',BloodFlow,BFFilt,Shear,SRFilt,Vel100',Diam100'];

% Write processed cyclic data to file
SumName1 = input('Enter filename to write cyclic data (must include file extension ".csv"): ', 's');
csvwrite(SumName1, AllData);

disp('Data Format: PSVTime,PSV,PDVTime,PDV,EDVTime,EDV,ISVTime,ISV,MBF1,MBF2,SysTime,DiasTime,OSI,WindowTime.');

% Write processed time series data to file
SumName2 = input('Enter filename to write time series data (must include file extension ".csv"): ','s');
csvwrite(SumName2,AllTimeSeries);

disp('Data Format: Time,BloodFlow,Filtered Blood Flow, Shear, Filtered Shear, Velocity, Diameter');

% Write analysis settings to file
SumName3 = input('Enter filename to write analysis settings (must include file extension ".csv"): ', 's');

% Store Epoch Selection Times
AnalysisSet = [PkHtSet', PkWdSet', RRDurKeep', PeakCount',DRoiTime,VelCal',TimeCal',DistCal',ZeroCal'];
csvwrite(SumName3,AnalysisSet);

disp('Analysis Settings Format: Peak Height Threshold, Peak Width Threshold, Cycle Duration, Peak Count, EpochEndTime, Velocity Calibration, Time Calibration, Distance Calibration, Zero Velocity Row Position');
