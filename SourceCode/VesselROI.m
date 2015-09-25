function [DROIx,DROIy,Center,Mask,theta] = VesselROI(I)
%VESSELROI define region of interest (ROI) around the vessel walls
%   VESSELROI allows the user to select a rectangular ROI around a vessl in
%   an ultrasound BMode image. 
%
%   [DROIx, DROIy, Center, Mask, theta] = VESSELROI(I) gets the image of
%   the ultrasound screen and returns user selected positions of the ROI x
%   and y coordinates, the center of the vessel, the mask width for the
%   vessel lumen, and the vessel wall angle. 
%
%   User selection is performed using the GINPUTC function. 
%
% Crystal Coolbaugh
% April 16, 2015
% Copyright 2015 Crystal Coolbaugh

VesselCheck = 0;

% Plot Image & Select ROI
[x,y] = ginputc(2,'Color','r','LineWidth',2);
DROIx = floor(x);
DROIy = floor(y);
close all;

while VesselCheck == 0
    
    % Select Center Point of Vessel
    hold on
    ICen = double(I(DROIy(1,1):DROIy(2,1),DROIx(1,1):DROIx(2,1)))';
    imagesc(ICen), axis image, colormap gray, axis off;title('Select center of vessel.');
    [X,Y] = ginputc(1,'Color','r','LineWidth',2);
    X=floor(X);
    Y=floor(Y);
    Center = [X Y]';
    close all;
    
    % Select Vessel Mask Width
    hold on
    imagesc(ICen), axis image, colormap gray, axis off;title('Select vessel mask width (2 points).');
    [X,Y] = ginputc(2,'Color','r','LineWidth',2);
    Mask = round(abs(X(2) - X(1)));
    plot(X,Y,'r','LineWidth',3);
    pause(0.5);
    close all;
    
    % Select Vessel Wall Angle
    hold on
    imagesc(ICen), axis image, colormap gray, axis off; title('Select 2 points on the vessel wall to define the angle');
    [X,Y] = ginputc(2,'Color','r','LineWidth',2);
    X=floor(X);
    Y=floor(Y);
    AngleOpp = X(1) - X(2);
    AngleAdj = abs(size(ICen,1) - 1);
    theta = [atand(AngleOpp/AngleAdj);1];
    plot(X,Y,'g','LineWidth',3);
    pause(0.5);
    close all;
    
    % Repeat Point Selection if Needed
    VesselTest = input('Were Vessel Points Selected Correctly? (Y/N)', 's');
    if lower(VesselTest) == 'y'
        disp('Advancing to Automated Diameter')
        close all;
        VesselCheck = 1;
    else
        disp('Reselect Vessel Points.')
        clear X Y;
    end
end