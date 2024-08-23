% Figure 1C 
% Written by: Josie D'Angelo 2024.

clear; close all;

% Setting path
current_path = pwd;
addpath(current_path);

if ismac
    cd([current_path, '/Data_Figure_1C']);
else
    cd([current_path, '\Data_Figure_1C']);
end

startDir = strcat( '', cd,'');
fileNames = dir('*.csv');

for i = 1: size(fileNames,1)
    
    trace_pixel = readtable(fileNames(i).name);
    
    % Setting default parameters
    expParameters = table(); expParameters.Subject = 1; expParameters.Gain = 0; expParameters.pathNum = 1;
    durationSec   = 1.5; % Duration in seconds
    frameRate     = size(trace_pixel,1)/durationSec; % Frames per second
    ppd_x         = 300; % Pixels per degree in the x dimension 
    ppd_y         = 300; % Pixels per degree in the y dimension 
    
    % Reformating for 'calculatingDandAlpha.m'
    trace_x_y = [trace_pixel.X_pixels, trace_pixel.Y_pixels];
    
    % Converting from pixels to degrees to arcminutes
    trace_arcmin(:,1) = trace_x_y(:,1)./ppd_x*60; % 300 pixels/deg. 1 deg = 60 arcmins;
    trace_arcmin(:,2) = trace_x_y(:,2)./ppd_y*60;
    
    % Defining variables
    overlapping     = 0; % Set to 0 because we used NONoverlapping intervals for all computations
    plotloggraph    = 0; % [optional]
    plotlinearGraph = 0; 
    plotIndiv       = 0; 
    
    % Computing alpha and diffusion constant [arcmin^2/s]
    [alpha, D] = calculatingDandAlpha(trace_arcmin, durationSec, expParameters, overlapping, plotloggraph, plotlinearGraph, plotIndiv);
    
    % Computing speed [arcmin/s]
    stepsDiff  = diff(trace_arcmin);
    hypotenuse = sqrt(sum(stepsDiff.^2,2));
    speed      = mean(hypotenuse).*frameRate; % Converting from frames to seconds
    
    % Plotting trace
    figure('color','w'); hold on;
    plot(trace_arcmin(:,1),trace_arcmin(:,2),'k-','linewidth',1.1); % Plotting trace
    plot(0,0,'+','color','r','markersize',10, 'linewidth',4); % Plotting starting position red cross
    plot([0 4],[-10 -10],'k-','linewidth',1.1); % Label for 4 arcmin bar
    axisLim = 20;
    set(gca,'xlim',[-axisLim axisLim], 'ylim',[-axisLim axisLim], 'Color', [1 1 1], 'FontSize', 18); axis square;
    sprintf('Fig: %1.f\nD: %1.f\nAlpha: %1.2f\nSpeed:  %1.f', i, D, alpha, speed)
end
cd ..; % Returning to ParadoxicalMisperception folder