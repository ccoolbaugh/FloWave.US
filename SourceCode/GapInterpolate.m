function [VelocityInterp] = GapInterpolate(Velocity, Time, TimeCon)
%GAPINTERPOLATE Linear interpolation scheme for linking individual frames
%
%   GAPINTERPOLATE uses a linear interpolation scheme to fill voids (200)
%   in an individual frame of data and to create a composite time series 
%   data set. 
%
%   [VelocityInterp] = GAPINTERPOLATE(Velocity, Time, TimeCon) gets the
%   velocity vector, the time vector, and the time pixel scaling factor and
%   returns an interpolated velocity vector. 
%

% Crystal Coolbaugh
% April 25, 2014
% Copyright 2014 Crystal Coolbaugh

clear GapIndex GoodIndex GapNum VelocityInterp

VelocityInterp = Velocity;
GapIndex = find(Velocity == 200);
GoodIndex = find(Velocity ~=200);
GapNum = length(GapIndex);


if isempty(GapIndex) ~= 1   % Interpolate if Gaps Exist
    for j = 1:GapNum
        
        % Find location of the gaps
        Low = find(GoodIndex < GapIndex(1));    % Data left of the gap
        High = find(GoodIndex > GapIndex(1));   % Data right of the gap
        
        % Data on both sides of gap
        if (isempty(High) == 0) && (isempty(Low) == 0)
            
            % Location of the points closest to the gap
            LowPos = GoodIndex(Low(end));
            HighPos = GoodIndex(High(1));
            
            % Find data for interpolation
            Y1 = VelocityInterp(LowPos);
            Y2 = VelocityInterp(HighPos);
            X1 = Time(LowPos);
            X2 = Time(HighPos);
            
            % Linear interpolation - find new point (X,Y)
            X = X1 + TimeCon;
            Y = Y1 + (Y2-Y1)*((X-X1)/(X2-X1));
            
            % Add new point into original vectors
            VelocityInterp(GapIndex(1),1) = Y;
            
            % Reset interpolation values
            clear X Y
            
            % Data on right side of gap
        elseif (isempty(High) == 0) && (isempty(Low) == 1)
            % Beginning of video - assign zero
            VelocityInterp(GapIndex(1),1) = 0;
            
            % Data on left side of gap
        elseif (isempty(High) == 1) && (isempty(Low) == 0)
            % Find data left of the gap
            LowPos = GoodIndex(Low(end));
            
            % Duplicate last good point
            VelocityInterp(GapIndex(1),1) = Velocity(LowPos);
            
            % No good data around the gap
        elseif (isempty(High) == 1) && (isempty(Low) == 0)
            % Assign zeros
            VelocityInterp(GapIndex(1),1) = 0;
        end
        
        % Find next gap
        GapIndex = find(VelocityInterp == 200);
        GoodIndex = find(VelocityInterp ~= 200);
    end
end