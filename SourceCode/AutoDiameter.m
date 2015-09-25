function DiameterPixel = AutoDiameter(I, Center, theta)
%AUTODIAMETER Measure vessel diameter using edge detection.
%   AUTODIAMETER uses the Canny edge detection function to detect the
%   vessel edges within a user defined region of interest (ROI). 
%
%   DiameterPixel = AUTODIAMETER(I, Center, theta) gets the image of the
%   user-defined ROI, the vessel center location, and the vessel wall angle
%   and returns the diameter in pixels for the specific image frame.  
%
%   This function is called by the program FloWave.US and the user-defined
%   ROI is created in the function VESSELROI.

% Crystal Coolbaugh
% April 16, 2015
% Copyright 2015 Crystal Coolbaugh

% Transpose Image
I = double(I)';

% Gaussian Filter
% Sigma = 1
% Kernel Size = 6*sigma
IFilt =  imgaussian(I,1);

% Edge Detect 
IEdge = edge(IFilt, 'canny');

NumRow = size(IEdge,1);
    
LowWall = zeros(NumRow,1);
UpWall = zeros(NumRow,1);
ImageDiameter = zeros(NumRow,1);

for i  = 1:NumRow
    loc = find(IEdge(i,:) == 1);
    if isempty(loc)
        ImageDiameter(i,1) = 200;
    else
        LowLoc = loc(find(loc < Center(1)));
        UpLoc = loc(find(loc > Center(1)));
        
        % Assign first pixel or identify gap
        if isempty(LowLoc) || isempty(UpLoc);
            ImageDiameter(i,1) = 200;
        else
            LowWall(i,1) = LowLoc(end);
            UpWall(i,1) = UpLoc(1);
            
            % Calculate Diameter - Theta Based on Vessel Wall Angle
            % Definition (User Input)
            ImageDiameter(i,1) = (UpWall(i,1) - LowWall(i,1))*cosd(theta);
        end
    end
end
    
ActivePixels = ImageDiameter(find(ImageDiameter ~= 200));

AvgD = mean(ActivePixels);

% Filter Erroneous Values
PixError = 0.10;    % Default = 10% difference from mean
if ~isnan(AvgD)
    index = 1;
    for i = 1:length(ActivePixels)
        if (ActivePixels(i,1)<(AvgD + PixError*AvgD)) && ...
                (ActivePixels(i,1)>(AvgD - PixError*AvgD))
            FilteredD(index,1) = ActivePixels(i,1);
            index = index+1;
        else
            FilteredD = 0;
        end
        
    end
    
    NonZeroFilterD = FilteredD(find(FilteredD ~= 0));
    
    if ~isempty(NonZeroFilterD)
        DiameterPixel = mean(NonZeroFilterD);  %Average diameter for each image
    else
        DiameterPixel = 200;
    end
else
    DiameterPixel = 200;
end

clear LowLoc UpLoc LowWall UpWall ImageDiameter FilteredD

