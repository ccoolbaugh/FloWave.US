function WriteFile(AllData,AllTimeSeries,AnalysisSet)

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
csvwrite(SumName3,AnalysisSet);

disp('Analysis Settings Format: Peak Height Threshold, Peak Width Threshold, Cycle Duration, Peak Count, EpochEndTime, Velocity Calibration, Time Calibration, Distance Calibration, Zero Velocity Row Position');


end

