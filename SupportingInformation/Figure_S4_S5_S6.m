% Figure S4, Figure S5, and Figure S6
% Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .

clear; close all; clc;

plotFigNum = 4; % 4--plotting S4;  5--plotting S5;  6--plotting S6

% Data folders for each figure
if plotFigNum == 4
    folder = '../SupportingInformation/Data_Figure_S4';
    subjectIDs = {'10003L'};
elseif plotFigNum == 5 
    folder = '../SupportingInformation/Data_Figure_S5';
    subjectIDs = {'10003L', '20114R'}; 
elseif plotFigNum == 6
    folder = '../SupportingInformation/Data_Figure_S6';
    subjectIDs = {'10003L'};
end

cd(folder);
 
for s = 1: length(subjectIDs)
    subjectID = num2str(subjectIDs{s});
    startDir = strcat( '', cd,'');
    fileNames = dir([startDir filesep subjectID '*.csv']);
    
    allGains = table();
    individGains = table();
    
    for i = 1: size(fileNames,1)
        
        % Table with averaged data
        if contains(fileNames(i).name, '_World_combined.csv') == 1
            cursubject = readtable(fileNames(i).name);
            allGains = [allGains; cursubject];
        end
        
        % Table with data from indivial trials
        if contains(fileNames(i).name, '_World_individTrial.csv') == 1
            block = readtable(fileNames(i).name);
            individGains = [individGains; block];
        end
    end
    
    % Defining background and Gain conditions
    numBackground = unique(allGains.Background); % "1" defines Background-present; "2" defines Background-absent
    numGains = unique(allGains.Gain);
    
    for b = 1: length(numBackground)
        
        f = figure('color','w'); hold on;
        [xLim] = formatGraphs(plotFigNum, b);

        % Color and shape labels
        color_plus = [0 0.35 0];
        color_minus = [0 0 0.6];
        if plotFigNum ~= 4
            color_third = [0.9290 0.6940 0.1250]; % Gain zero
        else
            color_third = [163/256 209/256 130/256]; % Orthogonal control exp 
        end
        shape = "o";

        for g = 1: length(numGains)
            % Selecting each Gain and background
            curGain_all = allGains(allGains.Gain == numGains(g) & allGains.Background == numBackground(b),:);
            curGain_indiv = individGains(individGains.Gain == numGains(g) & individGains.Background == numBackground(b),:);
            
            if numGains(g) == -1.5
                [~] = plotPoints(curGain_all, curGain_indiv, color_minus, shape, xLim);
            elseif numGains(g) == 1.5
                [~] = plotPoints(curGain_all, curGain_indiv, color_plus, shape, xLim);
            else 
                [~] = plotPoints(curGain_all, curGain_indiv, color_third, shape, xLim);
            end
        end
        
        % Plotting 1:1 line
        xmid = linspace(0,xLim,500);
        plot(xmid,xmid,'-', 'color', [0 0 0],'linewidth',1, 'HandleVisibility','off');
        
        xlabel('Diffusion constant for world motion, D_{WM} [arcmin^2/s]', 'FontSize', 14, 'fontname','helvetica');
        ylabel('Diffusion constant for perceived motion, D_{PM} [arcmin^2/s]', 'FontSize', 14, 'fontname','helvetica');
        set(gca, 'FontSize', 18, 'fontname','helvetica');
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});
        title(subjectID);
    end    
end

%% Functions
% Plotting individual subject points
function [xLim] = plotPoints(curGain_all, curGain_indiv, color, shape, xLim)
    sem_y = curGain_all.sem_D_Perceived;  sem_x= curGain_all.sem_D_World;
    if curGain_indiv.Gain(1) == -1.5 % Only Gain -1.5 hollow points
        errorbar(curGain_all.mean_D_World, curGain_all.mean_D_Perceived,sem_y,sem_y,sem_x, sem_x, shape,'Color', color, 'linewidth',1,'markersize',15);
    else
        errorbar(curGain_all.mean_D_World, curGain_all.mean_D_Perceived,sem_y,sem_y,sem_x, sem_x, shape,'Color', color, 'MarkerFaceColor', color, 'linewidth',1,'markersize',15);
    end
    
    [~]= makingQuivers(curGain_all.mean_Alpha_World, curGain_all.mean_D_World,curGain_all.mean_D_Perceived, xLim);

    for jp = 1: length(curGain_indiv.Gain)
        if curGain_indiv.Gain(1) == -1.5 % Only Gain -1.5 hollow points
            plot(curGain_indiv.D_World(jp), curGain_indiv.D_Perceived(jp), shape,'Color', color,'linewidth',1,'markersize',10);
        else
            plot(curGain_indiv.D_World(jp), curGain_indiv.D_Perceived(jp), shape,'Color', color, 'MarkerFaceColor', color,'linewidth',1,'markersize',10);
        end
        
        [~]= makingQuivers(curGain_indiv.Alpha_World(jp), curGain_indiv.D_World(jp),curGain_indiv.D_Perceived(jp), xLim);
    end
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
    qscale = quiver(xLim-xLim/2,xLim/4,(xLim * fulllength),0  ,0, 'Color', [1 0 0], 'LineWidth', 1.2, 'HandleVisibility','off'); 
    qscale.ShowArrowHead = 'off';
    plot(xLim-xLim/2+(xLim * fulllength),xLim/4,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');

    % Drawing arrow indicating anti/persistence
    if abs(alphaWorld-1) >0.02 && alphaWorld ~=0 
        q2 = quiver(worldMotion,perceivedmotion,quivLengthworld,0  ,0, 'Color', [1 0 0], 'LineWidth', 1.2, 'HandleVisibility','off');
        q2.ShowArrowHead = 'off';
        plot(worldMotion+quivLengthworld,perceivedmotion,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');
    end
end

% Formating axis
function [xLim] = formatGraphs(plotFigNum, b)
    grid on; axis square; box on;
    % Axes specific to subject due to idiosyncratic differences
    if plotFigNum == 4
        if b == 1
            xLim = 18;
            xticks(0:2:xLim); yticks(0:2:xLim);
            xlim([0 xLim]); ylim([0 xLim]);
        elseif b == 2
            xLim = 480;
            yLim = xLim/8;
            xticks(0:60:xLim); yticks(0:10:yLim);
            xlim([0 xLim]); ylim([0 yLim]);
        end
    elseif plotFigNum == 5
        xLim = 30; 
        xticks(0:5:xLim); yticks(0:5:xLim);
        xlim([0 xLim]); ylim([0 xLim]);
    elseif plotFigNum == 6
        xLim = 60;
        xticks(0:10:xLim); yticks(0:10/2:xLim/2);
        xlim([0 xLim]); ylim([0 xLim/2]);
    end  
end