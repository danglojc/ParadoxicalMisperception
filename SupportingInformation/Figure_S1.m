% Figure S1
% Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .

clear; close all;

% Setting path
current_path = pwd;
addpath(current_path);

if ismac
    cd([current_path, '/Data_Figures_2_3_6_S1']);
else
    cd([current_path, '\Data_Figures_2_3_6_S1']);
end
 
% Subject IDs
subjectIDs = {'10003L', '20114R', '20229L', '20234L', '20237R', '20256R'};
    
for s = 1:  length(subjectIDs)
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
        [xLim] = formatGraphs(b, subjectID);

        % Color and shape labels
        color_plus = [0 0.35 0];
        color_minus = [0 0 0.6];
        color_zero = [0.9290 0.6940 0.1250];
        shape = "o";

        for g = 1: length(numGains)
            % Selecting each Gain and Background
            curGain_all = allGains(allGains.Gain == numGains(g) & allGains.Background == numBackground(b),:);
            curGain_indiv = individGains(individGains.Gain == numGains(g) & individGains.Background == numBackground(b),:);
            
            if numGains(g) == -1.5
                [~] = plotPoints(curGain_all, curGain_indiv, color_minus, shape, xLim);
            elseif numGains(g) == 1.5
                [~] = plotPoints(curGain_all, curGain_indiv, color_plus, shape, xLim);
            elseif numGains(g) == 0
                [~] = plotPoints(curGain_all, curGain_indiv, color_zero, shape, xLim);
            end
        end
        
        % Plotting 1:1 line
        xmid = linspace(0,xLim,500);
        plot(xmid,xmid,'-', 'color', [0 0 0],'linewidth',1, 'HandleVisibility','off');
        
        xlabel('Diffusion constant for world motion, D_{WM} [arcmin^2/s]');
        ylabel('Diffusion constant for perceived motion, D_{PM} [arcmin^2/s]');
        set(gca, 'FontSize', 18, 'fontname','helvetica');
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'});
        title(subjectID);
    end   
end
cd ..; % Returning to ParadoxicalMisperception folder
%% Functions
% Plotting individual subject points
function [xLim] = plotPoints(curGain_all, curGain_indiv, color, shape, xLim)
    sem_y = curGain_all.sem_D_Perceived	;  sem_x= curGain_all.sem_D_World;
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
    qscale = quiver(xLim-xLim/2,xLim-xLim/2,(xLim * fulllength),0  ,0, 'Color', [1 0 0], 'LineWidth', 1.2, 'HandleVisibility','off'); 
    qscale.ShowArrowHead = 'off';
    plot(xLim-xLim/2+(xLim * fulllength),xLim-xLim/2,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');

    % Drawing arrow indicating anti/persistence
    if abs(alphaWorld-1) >0.02 && alphaWorld ~=0
        q2 = quiver(worldMotion,perceivedmotion,quivLengthworld,0  ,0, 'Color', [1 0 0], 'LineWidth', 1.2, 'HandleVisibility','off');
        q2.ShowArrowHead = 'off';
        plot(worldMotion+quivLengthworld,perceivedmotion,"diamond", 'Color', [1 0 0], 'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'HandleVisibility','off');
    end
end

% Formating axis
function [xLim] = formatGraphs(b, subjectID)
    numSub = str2double(subjectID(1:5));
    grid on; axis square; box on;
    if b == 1
        % Axes specific to subject due to idiosyncratic differences
        if numSub == 10003
            xLim = 10;
        elseif numSub == 20114
            xLim = 12;
        elseif numSub == 20229
            xLim = 16;
        elseif numSub == 20234
            xLim = 24;
        elseif numSub == 20237
            xLim = 20;
        elseif numSub == 20256
            xLim = 28;
        end
        
        if xLim <22
            xticks(0:2:xLim); yticks(0:2:xLim);
        else 
            xticks(0:4:xLim); yticks(0:4:xLim);
        end
        xlim([0 xLim]); ylim([0 xLim]);
    elseif b == 2
        % Axes specific to subject due to idiosyncratic differences
        if numSub == 10003 || numSub == 20237
            xLim = 140;
        elseif numSub == 20114
            xLim = 120;
        elseif numSub == 20229
            xLim = 50;
        elseif numSub == 20234 || numSub == 20256
            xLim = 70;
        end
        
        if numSub == 20234 || numSub == 20256
            yLim =  floor(xLim/2);
            xticks(0:10:xLim); yticks(0:10/2:yLim);
        else
            yLim =  (xLim/4);
            if xLim == 50
                xticks(0:10:xLim); yticks(0:10/4:yLim);
            else
                xticks(0:20:xLim); yticks(0:20/4:yLim);
            end 
        end
        xlim([0 xLim]);  ylim([0 yLim]);
    end
end