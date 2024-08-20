% Figure 2
% Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .

clear; close all; clc;

folder = '../Data_Figures_2_3_6_S1';
cd(folder);
 
% Table storing all subjects data
allGains = table();

startDir = strcat( '', cd,'');
fileNames = dir('*.csv');
for i = 1: size(fileNames,1)
    % Table with averaged data
    if contains(fileNames(i).name, '_World_combined.csv') == 1
        cursubject = readtable(fileNames(i).name);
        allGains = [allGains; cursubject];
    end  
end

% Color labels
color_plus  = [0 0.35 0];
color_minus = [0 0 0.6];
color_zero  = [0.9290 0.6940 0.1250]; 

% Defining variables
subjects      = unique(allGains.Subject);
numBackground = unique(allGains.Background); % "1" defines Background-present; "2" defines Background-absent
numGains      = unique(allGains.Gain);

for b = 1: length(numBackground)  

    f = figure('color','w'); hold on;
    [xLim] = formatGraphs(b);

    for s = 1:  length(subjects)
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
        
        for g = 1: length(numGains)
            % Selecting each Gain and Background
            curGain_all = allGains(allGains.Subject == subjects(s) & allGains.Gain == numGains(g) & allGains.Background == numBackground(b),:);
            
            if numGains(g) == -1.5
                [~] = plotSubjectPoints(curGain_all, color_minus, shape, xLim);
            elseif numGains(g) == 1.5
                [~] = plotSubjectPoints(curGain_all, color_plus, shape, xLim);
            elseif numGains(g) == 0
                [~] = plotSubjectPoints(curGain_all, color_zero, shape, xLim);
            end
        end       
    end
    
    % Plotting group averages
    [~] = plotAverages(b, allGains, subjects, numGains, color_plus, color_minus, color_zero);
    
    % Plotting 1:1 line
    xmid = linspace(0,xLim,500);
    plot(xmid,xmid,'-', 'color', [0 0 0],'linewidth',1, 'HandleVisibility','off');
    
    xlabel('Diffusion constant for world motion, D_{WM} [arcmin^2/s]', 'FontSize', 14, 'fontname','helvetica');
    ylabel('Diffusion constant for perceived motion, D_{PM} [arcmin^2/s]', 'FontSize', 14, 'fontname','helvetica');
    set(gca, 'FontSize', 18, 'fontname','helvetica');
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});
end
%% Functions
% Plotting group averages with standard error of the mean bars
function [color_zero] = plotAverages(b, allGains, subjects, numGains, color_plus, color_minus, color_zero)
    for g = 1:length(numGains)
        mean_x = mean(allGains.mean_D_World(allGains.Gain == numGains(g) & allGains.Background==b));
        sem_x = std(allGains.mean_D_World(allGains.Gain == numGains(g) & allGains.Background==b))/sqrt(length(subjects));
        mean_y = mean(allGains.mean_D_Perceived(allGains.Gain == numGains(g) & allGains.Background==b));
        sem_y = std(allGains.mean_D_Perceived(allGains.Gain == numGains(g) & allGains.Background==b))/sqrt(length(subjects));
        
        if numGains(g) == -1.5
            errorbar(mean_x, mean_y,sem_y,sem_y,sem_x, sem_x, "pentagram",'Color', color_minus, 'linewidth',1,'markersize',20);
        elseif numGains(g) == 1.5
            errorbar(mean_x, mean_y,sem_y,sem_y,sem_x, sem_x, "pentagram",'Color', color_plus, 'MarkerFaceColor', color_plus, 'linewidth',1,'markersize',20);    
        elseif numGains(g) == 0
            errorbar(mean_x, mean_y,sem_y,sem_y,sem_x, sem_x, "pentagram",'Color', color_zero, 'MarkerFaceColor', color_zero, 'linewidth',1,'markersize',20);
        end
    end
end

% Plotting individual subject points
function [xLim] = plotSubjectPoints(curGain_all, color, shape, xLim)  
    if curGain_all.Gain(1) == -1.5 % Only Gain -1.5 hollow points
        plot(curGain_all.mean_D_World, curGain_all.mean_D_Perceived, shape,'Color', color, 'linewidth',1,'markersize',9);
    else
        plot(curGain_all.mean_D_World, curGain_all.mean_D_Perceived, shape,'Color', color, 'MarkerFaceColor', color, 'linewidth',1,'markersize',9);
    end
    
    [~]= makingQuivers(curGain_all.mean_Alpha_World, curGain_all.mean_D_World,curGain_all.mean_D_Perceived, xLim);
end

% Adding red alpha quivers to indicate anti/persistence
function [alphaWorld]= makingQuivers(alphaWorld, worldMotion,perceivedmotion, xLim)
    
    fulllength = 20/140; % So that the length will be the same regardless of axis

    % Alpha World motion
    if alphaWorld ~= 0
        quivLengthworld = (alphaWorld-1)* xLim * fulllength;
    else
        quivLengthworld = 0;
    end

    % Plotting a quiver length that is equivalent to alpha = 2 or alpha = 0, used for the legend
    qscale = quiver(xLim-xLim+2,xLim-xLim+2,(xLim * fulllength),0  ,0, 'Color', [1 0 0], 'LineWidth', 1.2, 'HandleVisibility','off');
    qscale.ShowArrowHead = 'off';
    plot(xLim-xLim+2+(xLim * fulllength),xLim-xLim+2,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');

    % Drawing arrow indicating anti/persistence
    if abs(alphaWorld-1) >0.02 && alphaWorld ~=0
        q2 = quiver(worldMotion,perceivedmotion,quivLengthworld,0  ,0, 'Color', [1 0 0], 'LineWidth', 1.2, 'HandleVisibility','off');
        q2.ShowArrowHead = 'off';
        plot(worldMotion+quivLengthworld,perceivedmotion,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');
    end

end

% Formating axis
function [xLim] = formatGraphs(b)
    grid on; axis square; box on;
    if b == 1
        xLim = 20;
        xticks(0:2:xLim); yticks(0:2:xLim);
        xlim([0 xLim]); ylim([0 xLim/2]);
    elseif b == 2
        xLim = 100;
        xticks(0:20:xLim); yticks(0:20/4:xLim/4);
        xlim([0 xLim]); ylim([0 xLim/4]);
    end
end