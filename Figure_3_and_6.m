% Figure 3 and Figure 6
% Written by: Josie D'Angelo 2024.
clear; close all; clc;

plotFigNum = 3; % 3--plotting Figure 3;   6--plotting Figure 6

% Setting path
current_path = pwd;
addpath(current_path);

if ismac
    cd([current_path, '/Data_Figures_2_3_6_S1']);
else
    cd([current_path, '\Data_Figures_2_3_6_S1']);
end

startDir = strcat( '', cd,''); 
fileNames = dir([startDir filesep '*.csv']); 

all_subects_allgain = table();

for i = 1: size(fileNames,1) 
    % Table with averaged data
    if contains(fileNames(i).name, '_World_combined.csv') == 1 
        cursubject = readtable(fileNames(i).name);
        all_subects_allgain = [all_subects_allgain; cursubject];
    end  
end

% Defining variables
numBackground = unique(all_subects_allgain.Background); % "1" defines Background-present; "2" defines Background-absent
numGains      = unique(all_subects_allgain.Gain);
subjects      = unique(all_subects_allgain.Subject);

% Ratios: Row 1 stores Gain -1.5 ratios and Row 2 stores Gain +1.5 ratios
ratioBAbsentX  = nan(2,length(subjects));
ratioBPresentY = nan(2,length(subjects));  

% Alphas: Row 1 stores Gain -1.5 alphas and Row 2 stores Gain +1.5 alphas
alphaWorldBAbsentX  = nan(2,length(subjects)); 
alphaWorldBPresentY = nan(2,length(subjects));

% Quiver lengths: Row 1 stores Gain -1.5 quiver lengths and Row 2 stores Gain +1.5 quiver lengths
quivLengthX = nan(2,length(subjects)); 
quivLengthY = nan(2,length(subjects));

% Setting up Figure 3
figure('color','w'); hold on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.8, 0.8]);

element = 1;
axisLim = 1.5;
labels = {};

color_plus = [0 0.35 0];
color_minus = [0 0 0.6];

for g = 1: length(numGains)

    for s = 1 : length(subjects)
        
        % Shape labels
        if s == 1
            shape = "o";
        elseif s == 2
            shape = "square";
        elseif s == 3
            shape = "^";
        elseif s == 4
            shape = "v";
        elseif s == 5
            shape = "hexagram";
        elseif s == 6
            shape = "diamond";
        end
        
        data_of_interest = all_subects_allgain(all_subects_allgain.Gain == numGains(g) & all_subects_allgain.Subject == subjects(s),:);
        
        if numGains(g) == -1.5
            row_gain = 1;
            [ratioBPresentY, ratioBAbsentX,alphaWorldBPresentY, alphaWorldBAbsentX,quivLengthX, quivLengthY, element, labels]= Figure_3_plotRatios(data_of_interest, ratioBPresentY, ratioBAbsentX, alphaWorldBPresentY, alphaWorldBAbsentX, axisLim, quivLengthX, quivLengthY, s, color_minus, shape, row_gain, element, labels);
            
        elseif numGains(g) == 1.5
            row_gain = 2;
            [ratioBPresentY, ratioBAbsentX,alphaWorldBPresentY, alphaWorldBAbsentX,quivLengthX, quivLengthY, element, labels]= Figure_3_plotRatios(data_of_interest, ratioBPresentY, ratioBAbsentX, alphaWorldBPresentY, alphaWorldBAbsentX, axisLim, quivLengthX, quivLengthY, s, color_plus, shape, row_gain, element, labels);
            
        end
    end
end

% Plotting Figure 6
if plotFigNum == 6
    [~] = Figure_6_plotRatioVsAlpha(subjects, ratioBAbsentX, alphaWorldBAbsentX, ratioBPresentY, alphaWorldBPresentY, color_plus, color_minus);
    cd ..; % Returning to ParadoxicalMisperception folder
    return;
end

% Plotting Figure 3
% Plotting yellow lines connecting the points
for s = 1: length(subjects)
    plot([ratioBAbsentX(2, s) ratioBAbsentX(1, s)],[ratioBPresentY(2, s) ratioBPresentY(1, s)], 'Color', [0.9290 0.8940 0.1250], 'linewidth',1, 'HandleVisibility','off');
end

% Adding red alpha quivers to indicate anti/persistence
[~] = makeQuivers(alphaWorldBAbsentX, alphaWorldBPresentY, ratioBPresentY, ratioBAbsentX, quivLengthX, quivLengthY);

% Plotting 1:1 line
xmid = linspace(0,axisLim,500);
plot(xmid,xmid,'-', 'color', [0 0 0],'linewidth',1, 'HandleVisibility','off') 

grid on; axis square; box on;
xlabel('[D_{pm} (arcmin^2/s) / D_{wm} (arcmin^2/s)]', 'FontSize', 14);  
ylabel('[D_{pm} (arcmin^2/s) / D_{wm} (arcmin^2/s)]', 'FontSize', 14);
xlim([0 (axisLim)]); ylim([0 (1.5)]);
xticks(0:0.25:axisLim); yticks(0:0.25:axisLim);
set(gca, 'FontSize', 20, 'fontname','helvetica')
    
h = get(gca,'Children');
legend(labels,'FontSize', 22);
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

cd ..; % Returning to ParadoxicalMisperception folder
%% Functions for Figure 3
% Plotting background-present ratios (y-axis) vs background-absent ratios (x-axis), and saving quiver lengths to indicate anti/persistence
function [ratioBPresentY, ratioBAbsentX,alphaWorldBPresentY, alphaWorldBAbsentX,quivLengthX, quivLengthY, element, labels]= Figure_3_plotRatios(data_of_interest, ratioBPresentY, ratioBAbsentX,alphaWorldBPresentY, alphaWorldBAbsentX, axisLim, quivLengthX, quivLengthY, s, color, shape,row_gain, element, labels)

    % Y axis: Background-present
    curBackground                    = data_of_interest(data_of_interest.Background == 1,:); % Background-present defined by "data_of_interest.Background == 1" 
    ratioBPresentY(row_gain, s)      = curBackground.mean_D_Perceived/curBackground.mean_D_World;
    alphaWorldBPresentY(row_gain, s) = curBackground.mean_Alpha_World;

    labels{element} = sprintf('%1.f:  %1.1f',data_of_interest.Subject(1), data_of_interest.Gain(1));
    element = element+1;

    % X axis: Background-absent
    curBackground                   = data_of_interest(data_of_interest.Background == 2,:); % Background-absent defined by "data_of_interest.Background == 2" 
    ratioBAbsentX(row_gain, s)      = curBackground.mean_D_Perceived/curBackground.mean_D_World;
    alphaWorldBAbsentX(row_gain, s) = curBackground.mean_Alpha_World;
    
    if color == [0 0.35 0]  
        plot(ratioBAbsentX(row_gain, s), ratioBPresentY(row_gain, s), shape,'Color', color, 'MarkerFaceColor', color,'linewidth',1,'markersize',15); hold on;
    else
        plot(ratioBAbsentX(row_gain, s), ratioBPresentY(row_gain, s), shape,'Color', color,'linewidth',1,'markersize',15); hold on;
    end

    fulllength = 15/140; % So that the length will be the same regardless of axis
    
    % X axis: Background-absent
    quivLengthX(row_gain, s) = (alphaWorldBAbsentX(row_gain, s)-1)* axisLim * fulllength;
    
    % Y axis: Background-present
    quivLengthY(row_gain, s) = (alphaWorldBPresentY(row_gain, s)-1)* axisLim * fulllength;

    % Plotting a quiver length that is equivalent to alpha = 2 or alpha = 0, used for the legend
    qscale = quiver(axisLim-axisLim/2,axisLim-axisLim/2,(axisLim * fulllength),0  ,0, 'Color', [1 0 0], 'LineWidth', 1.5, 'HandleVisibility','off'); 
    qscale.ShowArrowHead = 'off';
    plot(axisLim-axisLim/2+(axisLim * fulllength),axisLim-axisLim/2,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off')

end

% Adding red alpha quivers to indicate anti/persistence
function [ratioBAbsentX] = makeQuivers(alphaWorldBAbsentX, alphaWorldBPresentY, ratioBPresentY, ratioBAbsentX, quivLengthX, quivLengthY)
    for s = 1: size(quivLengthX,2) % For each subject
        for g = 1: size(quivLengthX,1) % For Gain -1.5 (g=1) and Gain +1.5 (g=2)
            
            % X axis: Background-absent
            if abs(alphaWorldBAbsentX(g,s)-1) > 0.02
                q1 = quiver(ratioBAbsentX(g,s), ratioBPresentY(g,s),0,quivLengthY(g,s),  0, 'Color', [1 0 0], 'LineWidth', 1.5, 'HandleVisibility','off'); 
                q1.ShowArrowHead = 'off';
                plot(ratioBAbsentX(g,s),ratioBPresentY(g,s)+quivLengthY(g,s),"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');
            end
            % Y axis: Background-present
            if abs(alphaWorldBPresentY-1) > 0.02
                q2 = quiver(ratioBAbsentX(g,s), ratioBPresentY(g,s),quivLengthX(g,s),0  ,0, 'Color', [1 0 0], 'LineWidth', 1.5, 'HandleVisibility','off'); 
                q2.ShowArrowHead = 'off';
                plot(ratioBAbsentX(g,s)+quivLengthX(g,s),ratioBPresentY(g,s),"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off')
            end
        end
    end
end

%% Function for Figure 6 
% Plotting two different conditions of ratios (y-axis) vs alpha for world motion (x-axis), and adding model
function [color_plus] = Figure_6_plotRatioVsAlpha(subjects, ratioBAbsentX, alphaWorldBAbsentX, ratioBPresentY, alphaWorldBPresentY, color_plus, color_minus)
    close all;

    element = 1;
    labels = [];

    % Setting up Figure 6
    figure('color','w'); hold on;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.8, 0.8]);

    % Plotting Gain -1.5 background-absent ratios
    [labels, element] = plotRatiovsAlpha(alphaWorldBAbsentX(1,:),ratioBAbsentX(1,:),color_minus, element, subjects, labels);
    
    % Plotting Gain +1.5 background-present ratios
    [labels, element] = plotRatiovsAlpha(alphaWorldBPresentY(2,:),ratioBPresentY(2,:),color_plus, element, subjects, labels); 
    
    % Modelling ratios vs Alpha relationship
    alphaXvals= linspace(0.9,1.8,1000);
    model = (60).^(1-alphaXvals); 
    plot(alphaXvals, model, 'Color', [0 0 0], 'Linewidth', 2); 
    labels{element} = strcat('D_{PM}/D_{WM} = 60^{1-α}');
  
    grid on; axis square; box on;
    xlim([0.9 1.8]); xticks(0.9:0.1:1.8);
    ylim([0 1.5]); yticks(0:0.25:1.5);
    xlabel('\alpha for world motion, \alpha_{WM}'); ylabel('[D_{PM} (arcmin^2/s) / D_{WM} (arcmin^2/s)]');
    set(gca, 'FontSize', 20, 'fontname','helvetica');
    
    h = get(gca,'Children');
    legend(labels,'FontSize', 22);
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});

    function [labels, element] = plotRatiovsAlpha(alphas, ratios, color, element, subjects, labels)
        for s = 1: length(subjects)
            if s == 1
                shape = "o";
            elseif s == 2
                shape = "square";
            elseif s == 3
                shape = "^";
            elseif s == 4
                shape = "v";
            elseif s == 5
                shape = "hexagram";
            elseif s == 6
                shape = "diamond";
            end

            if color == [0 0.35 0]
                plot(alphas(s), ratios(s), shape,'Color', color, 'MarkerFaceColor', color,'linewidth',1,'markersize',15); hold on;
                labels{element} = sprintf('%1.f:  %1.1f',subjects(s), 1.5);
            else % Only Gain -1.5 hollow points
                plot(alphas(s), ratios(s), shape,'Color', color,'linewidth',1,'markersize',15); hold on;
                labels{element} = sprintf('%1.f:  %1.1f',subjects(s), -1.5);
            end
            element =element+1;
        end

    end
end