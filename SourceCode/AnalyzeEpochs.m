function [AllPSVTime,AllPSV,AllPDVTime,AllPDV,AllEDVTime,AllEDV,AllMBF1,AllData,AllTimeSeries,AnalysisSet] = AnalyzeEpochs(EpochTime,Time100,DRoiTime,BFFilt,SRFilt,Vel100,Diam100,FsNew,VelCon,TimeCon,DistCon,BaseY1,BloodFlow,Shear)
%%% Cardiac Cycle Analysis
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
            WinBF = BF(Left(1):Right(end));
            WinSR = SR(Left(1):Right(end));
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
            
            % Quality Check - Three Cycles
            % If points are poorly selected, a new frequency should be
            % used and analysis repeated (BETA). 
            if (numel(DblPeaks) >= 3)
                if i <= 3
                    if PDVLoc ~= 0
                        KeyPts = [ISVLoc, PSVLoc, PDVLoc,EDVLoc];
                    else
                        KeyPts = [ISVLoc,PSVLoc,EDVLoc];
                    end
                    %{
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
                    %}
                elseif i == numel(DblPeaks)                
                        AnalyzeSuccess = AnalyzeSuccess + 1;
                end
            else
                AnalyzeSuccess = AnalyzeSuccess+1;
            end
            
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

%%% Create Composite Data Set for All Epochs

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

% Aggregate Cardiac Cycle Data  - Format for Validation Study Only
AllData = [AllPSVTime',AllPSV',AllPDVTime',AllPDV',AllEDVTime',AllEDV', ...
    AllISVTime',AllISV',AllMBF1',AllMBF2',...
    AllSysTime',AllDiasTime',AllOSI',AllWindow'];
 
% Aggregate Time Series Data
AllTimeSeries = [Time100',BloodFlow,BFFilt,Shear,SRFilt,Vel100',Diam100'];

% Store Epoch Selection Times
AnalysisSet = [PkHtSet', PkWdSet', RRDurKeep', PeakCount',DRoiTime,VelCal',TimeCal',DistCal',ZeroCal'];



end

