function [TimeFinal,VelInterp] = CompositePWData(TimeCon,VelAll,MeanVelAll,EnvelopeVelAll,WtMeanVelAll,VelType)
%%% Create Composite Pulse-Wave Data
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


end

