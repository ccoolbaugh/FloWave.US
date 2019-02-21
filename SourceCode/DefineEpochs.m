function [EpochTime,DRoiTime,DFilt,DTime] = DefineEpochs(TimeFinal,VelInterp,TimeCon,VidSamRate,nFrames,USObj,DistCon)
%%% Define Data Epochs (VesselROI, AutoDiameter, imgaussian functions)
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


end

