function [VelType,VelAll,MeanVelAll,WtMeanVelAll,EnvelopeVelAll] = PulseExtraction(StartFrame,FramesPerPeriod,nFrames,USObj,BaseY1,VelCon,PlatSetFile)
%%% Pulse Wave Data Extraction 
% For all frames in the video, this code a) defines an ROI for the pulse
% wave data, b) extracts the time average mean data (calculated teal line),
% c) extracts the max doppler envelope (TAMax), and e) scales 
% the velocity and time data from pixels to the appropriate units.

% NOTE: If time average mean data is not displayed in real-time, the mean 
% of the doppler envelope may be used as an alternative. 

load(PlatSetFile)

% Initialize Variables that may not be necessary so the function can export
% all required values
VelAll = [];
MeanVelAll = [];
WtMeanVelAll = [];
EnvelopeVelAll = [];

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


end

