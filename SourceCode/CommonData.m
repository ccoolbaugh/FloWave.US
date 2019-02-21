function [FsNew,Time100,Vel100,Diam100,BloodFlow,Shear,BFFilt,SRFilt] = CommonData(VidSamRate,TimeCon,DTime,DFilt,TimeFinal,VelInterp)
%%% Create Common Time Series - Velocity & Diameter / Calculate Blood Flow and Shear Rate
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
%%% Calculate Blood Flow and Shear Rate

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


end

