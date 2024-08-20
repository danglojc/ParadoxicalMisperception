function [alpha, D_log, msd_each_path] = calculatingDandAlpha(paths_x_y, durationSec, expParameters, overlapping, plotloggraph, plotlinearGraph, plotIndiv)
    
    % Computes the log10(mean square displacement) vs log10(time interval) for all paths, output: the diffusion constant (D) and alpha 
    %    -- Equation: MSD = 2dD(dT)^alpha 
    %           -- MSD: mean square displacement,  d: dimension,  D: diffusion constant,  dT: time interval,  alpha: scaling exponent 
    %    -- Please note: Diffusion constant and alpha are computed from the log-log graph 
    %    -- Units: Input: 'paths_x_y' [arcmin]. Output: 'D_log' [arcmin^2/s], alpha [au]
    %     
    %   -----------------------------------
    %   Input
    %   ----------------------------------- 
    %   paths_x_y               : 3D array: positions in x--column1 [arcmin] and y--column2 [arcmin], multiple paths in 3rd axis 
    %   durationSec             : stimulus duration in seconds (ie 0.75 seconds)
    %   expParameters           : table with experiment parameters, each row corresponds to each path 
    %                              (default: expParameters.Subject--subject ID, expParameters.Gain--retina-contingent condition, expParameters.pathNum--indexs the path number) 
    %   overlapping             : defines whether to compute the mean square displacement across overlapping or NONoverlapping intervals
    %                              (1--overlapping intervals, 0--NONoverlapping intervals)
    %   plotloggraph            : plots average log10(mean square displacement) vs log10(time interval) across ALL paths from paths_x_y 
    %                              (1--plot figure, 0--don't plot)   
    %   plotlinearGraph         : plots average (mean square displacement) vs (time interval) across ALL paths from paths_x_y
    %                              (1--plot figure, 0--don't plot) 
    %   plotIndiv               : plots log10(mean square displacement) vs log10(time interval) for each individual path
    %                              (1--plot figure(s), 0--don't plot) 

    %   -----------------------------------
    %   Output
    %   -----------------------------------
    %   alpha                   : the scaling exponent computed across ALL paths from paths_x_y
    %   D_log                   : the diffusion constant computed across ALL paths from paths_x_y
    %   msd_each_path           : lists the D and alpha for each individual path 
    %                             (default: expParameters.Subject--subject ID, expParameters.Gain--retina-contingent condition, expParameters.pathNum--indexs the path number, expParameters.DiffusionConstant_indiv--D for individual path, expParameters.AlphaIndiv--alpha for individual path) 
    %
    % Associated with "A paradoxical misperception of relative motion" Josephine C. D’Angelo, Pavan Tiruveedhula, Raymond J. Weber, David W. Arathorn, Austin Roorda
    % Used to compute the diffusion constants and alphas for world motion, eye motion, and perceived motion
    % Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .


    % Setting up for the mean square displacement calculation
    numFrames     = size(paths_x_y,1);  % Number of frames
    frameRate     = numFrames/durationSec; % Frame rate of the system (ie 60 Hz)
    totNumPaths   = size(paths_x_y,3);  % Best practice to have multiple paths (I usually aim for ~20 paths)
    max_numDeltaT = floor(numFrames/4); % Maximum number of time intervals, recommended by Saxton et al.,1997--https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1184368/pdf/biophysj00037-0260.pdf
    all_msd       = nan(max_numDeltaT,4,totNumPaths);  % To track the (mean square displacements) per path
    all_logmsds   = nan(max_numDeltaT,4,totNumPaths); % To track the log10(mean square displacements) per path
 
    % Part 1: Calcuating (mean square displacement) and log10(mean square displacement) from x,y values of each path-------------------------------------------------------------------------------------- 
    for i = 1:totNumPaths
        data = paths_x_y(:,:, i);

        numberOfdeltaT = floor(numFrames/4); % time intervals up to 1/4 of number of total points (Saxton et al., 1997)
 
        msd     = nan(numberOfdeltaT,4); % store [mean, std, n, timelag[frames]]
        logmsds = nan(numberOfdeltaT,4); % store [log10(mean), log10(std), n, log10(timelag[frames])]
        
        % Calculate msd for all deltaT's, modified from https://stackoverflow.com/questions/7489048/calculating-mean-squared-displacement-msd-with-matlab
        for dt = 1:numberOfdeltaT
 
            if overlapping == 1
                deltaCoords = data(1+dt:end,1:2) - data(1:end-dt,1:2); %overlapping time intervals
            elseif overlapping == 0
                deltaCoords = data(1+dt:dt:end,1:2) - data(1:dt:end-dt,1:2); %NONoverlapping time intervals
            end
            
            squaredDisplacement = sum(deltaCoords.^2,2); % Units are arcmin^2 here
            
            % Linear plots
            msd(dt,1) = mean(squaredDisplacement); % average
            msd(dt,2) = std(squaredDisplacement);  % standard deviation
            msd(dt,3) = length(squaredDisplacement); % number of squared displacements 
            msd(dt,4) = dt; % time intervals (units in frames)
          
            % Log-log plots
            logmsds(dt,1) = log10(mean(squaredDisplacement)); % log10 average
            logmsds(dt,2) = log10(std(squaredDisplacement));  % log10 standard deviation
            logmsds(dt,3) = length(squaredDisplacement); % number of squared displacements 
            logmsds(dt,4) = log10(dt); % log10 time intervals 
        end
 
        % Making figures for individual paths
        if plotIndiv == 1
            [~] = plotIndivGraphs_log(logmsds, durationSec, numFrames, overlapping, i, frameRate);
        end
        
        all_msd(1:numberOfdeltaT,:,i)     = msd;
        all_logmsds(1:numberOfdeltaT,:,i) = logmsds;
    end
    
    % Calculating the average across all paths
    islinear = 1;
    [meanMSD] =  calcMeanAcrossAllPaths(max_numDeltaT, all_msd, islinear);

    islinear = 0;
    [meanlogMSD] =  calcMeanAcrossAllPaths(max_numDeltaT, all_msd, islinear);
    
    % Part 2: Plottting (mean square displacement) vs (time interval) across ALL paths-----------------------------------------------------------------------
    if plotlinearGraph == 1
        figure('color','w'); hold on;
        errorbar(meanMSD(:,4), meanMSD(:,1),meanMSD(:,2), 'o');
        DC = polyfit(meanMSD(:,4),meanMSD(:,1),1); 
        
        newX = linspace(0,65,10000);
        yline = DC(1)*newX+DC(2);
        [~, indexforclosestFrameRt] = min(abs(newX - frameRate)); % To convert from frames to seconds
        dc_linear = yline(indexforclosestFrameRt); 
        dc_linear = dc_linear/4; % MSD = 2dD(dT)^alpha, the final step is to divide by 2d. JD added 9/27/23
 
        onlydataptsYline = DC(1)*meanMSD(:,4)+DC(2);
        plot(meanMSD(:,4),onlydataptsYline,'r-');
        if overlapping == 1 
            title([sprintf('Gain %1.1f',unique(expParameters.Gain)), '. AVG (', num2str(totNumPaths), ' paths), overlapping Intervals']);
        elseif overlapping == 0
            title([sprintf('Gain %1.1f',unique(expParameters.Gain)), '. AVG (' , num2str(totNumPaths), ' paths), NONoverlappping Intervals']);
        end
        text(0.05, 0.25, sprintf('D: %.4f arcmin^2/s', dc_linear), 'FontSize', 16);
        xlabel('Time interval [milliseconds]');
        ylabel('Mean square displacement [arcmin^2]');
        set(gca, 'FontSize', 16);
        [~] = formatLineargraph(durationSec,numFrames);
    end
    
    % Part 3: Plotting log10(mean square displacement) vs log10(time interval) across ALL paths to get alpha and diffusion constant-----------------------------
    [D_log, alpha, logM] = interpXvalsgetDAlpha(meanlogMSD, plotloggraph, frameRate);
    
    if plotloggraph == 1 
        [posDlog, posAlphalog] = formatLoggraph(durationSec, numFrames);
        
        if exist('posAlphalog', 'var') && exist('posDlog', 'var')
            text(posAlphalog(1), posAlphalog(2), ['\alpha: ', sprintf('%.4f', alpha)], 'FontSize', 16);
            text(posDlog(1), posDlog(2), sprintf('D: %.4f arcmin^2/s', D_log), 'FontSize', 16);
        else
            text(0.1, 1.2, ['\alpha: ', sprintf('%.4f', alpha)], 'FontSize', 16);
            text(0.1, 1.4, sprintf('D: %.4f arcmin^2/s', D_log), 'FontSize', 16);
        end
        
        if overlapping == 1
            title([sprintf('Gain %1.1f',unique(expParameters.Gain)), '. AVG (', num2str(totNumPaths), ' paths), overlapping Intervals, LOG-LOG']);
        elseif overlapping == 0
            title([sprintf('Gain %1.1f',unique(expParameters.Gain)), '. AVG (' , num2str(totNumPaths), ' paths), NONoverlappping Intervals, LOG-LOG']);
        end
        xlabel('log_1_0 (Time interval) [milliseconds]');
        ylabel('log_1_0(MSD) [arcmin^2]');
        set(gca, 'FontSize', 16)
    end

    % Part 4: Calculating diffusion constant & alpha for each individual path JD 12/18/22----------------------------------------------------------------------
    msd_each_path = table();
    for t = 1: size(all_msd,3)
        msd = all_logmsds(:,:,t);
        
        % Saving the diffusion constant computed from the log-log graph
        [DC_of_interest, alpha_indiv] = interpXvalsgetDAlpha(msd, 0, frameRate);

        msd_each_path.Subject(t) = expParameters.Subject(t);
        msd_each_path.Gain(t)    = expParameters.Gain(t);
        msd_each_path.pathNum(t) = expParameters.pathNum(t);
        msd_each_path.DiffusionConstant_indiv(t) = DC_of_interest;
        msd_each_path.AlphaIndiv(t) = alpha_indiv;
    end
end
 
%% All Functions

% Plotting each individual path's log-log graph
function [DC_of_interest] =  plotIndivGraphs_log(logmsds, durationSec, numFrames, overlapping, i, frameRate)

    [DC_of_interest, alpha_indiv] = interpXvalsgetDAlpha(logmsds, 1, frameRate);
 
    [posDlog, posAlphalog] = formatLoggraph(durationSec, numFrames);  
    if exist('posAlphalog', 'var') && exist('posDlog', 'var')
        text(posAlphalog(1), posAlphalog(2), ['\alpha: ', sprintf('%.4f', alpha_indiv)], 'FontSize', 14)
        text(posDlog(1), posDlog(2), sprintf('D: %.4f arcmin^2/s', DC_of_interest), 'FontSize', 14)
    else
        text(0.1, 1.2, ['\alpha: ', sprintf('%.4f', alpha_indiv)], 'FontSize', 14)
        text(0.1, 1.4, sprintf('D: %.4f arcmin^2/s', DC_of_interest), 'FontSize', 14)
    end

    xlabel('log_1_0 (Time interval) [milliseconds]');
    ylabel('log_1_0(MSD) [arcmin^2]');
    
    if overlapping == 1
        title(sprintf('Path: %1.f log-log (overlapping intervals)', i));
    elseif overlapping == 0
        title(sprintf('Path: %1.f log-log (NONoverlap intervals)', i));
    end
 
end
 
% Computing the Average across all paths 
    % islinear == 1 means original msd 
    % islinear ~= 1 means convert to logspace: Compute the log of the mean (not the mean of the logs)
function [meanMSD] =  calcMeanAcrossAllPaths(max_numberOfdeltaT, all_paths, islinear)
    for row = 1:max_numberOfdeltaT
        for column = 1:4
            if islinear == 1
                if column ~= 2 % Since column 2 is stdev column (and don't want to take the mean of the stdevs)
                    meanMSD(row,column) = mean(all_paths(row, column, :));
                else
                    meanMSD(row,column) = std(all_paths(row, 1, :));
                end
            else % log data instead
                if column ~= 2
                    meanMSD(row,column) = log10(mean(all_paths(row, column, :)));
                else
                    meanMSD(row,column) = log10(std(all_paths(row, 1, :)));
                end
            end
        end
    end
end
 
% Interpolating the log-log graph to add weight to the lower x-values
function [D_log, alpha, logM] = interpXvalsgetDAlpha(meanlogMSD, plotgraph, frameRate)   
    
    % Plotting the log-log graph
    convertSD_up   = log10(10.^meanlogMSD(:,1)+ 10.^meanlogMSD(:,2)) - meanlogMSD(:,1);    %converting error bars since now in log-log space
    convertSD_down = meanlogMSD(:,1) - (log10(10.^meanlogMSD(:,1) - 10.^meanlogMSD(:,2))); %converting error bars since now in log-log space
    if plotgraph == 1
        figure('color','w'); hold on;
        errorbar(meanlogMSD(:,4), meanlogMSD(:,1),convertSD_down,convertSD_up, 'o');
    end

    % Interpolating lower x-values
    distanceX = diff(meanlogMSD(:,4)); % Interpolation weights are dependent on the log10(time intervals) spacing
    totnum = sum(distanceX);
    weight = distanceX./totnum; % Weights sum to 1
    
    all_intX = [];
    all_intY = []; 
    for pt = 1: length(meanlogMSD(:,4))-1
        intX = linspace(meanlogMSD(pt,4),meanlogMSD(pt+1,4), 1000*weight(pt)); % I chose 1000 ad hoc
        intY = interp1(meanlogMSD(pt:pt+1,4),meanlogMSD(pt:pt+1,1),intX,'linear'); % One-dimensional data interpolation 
        
        if isempty(intX) % If the weights to the right are zero, then you wouldn't interpolate and just save the original value 
            intX = meanlogMSD(pt,4);
            intY = meanlogMSD(pt,1);
        end
        
        all_intX = [all_intX, intX];
        all_intY = [all_intY, intY]; 
        if plotgraph == 1
            plot(intX,intY,'r.');  hold on;
        end
    end

    logM = polyfit(all_intX,all_intY,1);  % Getting line of best fit across interp log10(mean square displacement) vs log10(time interval)
    logmYline = logM(1)*all_intX+logM(2); % (polyval)
    
    if plotgraph == 1
        plot(all_intX,logmYline,'g-', 'LineWidth', 2)
    end
    
    % The time intervals are in units of [frames].
    % However, we want the diffusion constant in units of [arcmin^2/s]. 
    % The framerate of the AOSLO is 60 Hz, so we evaluted at the framerate which corresponds to 1 second
    newX = linspace(0,2,10000); % We determined log10(mean square displacement) at 60 frames 
    newYline = logM(1)*newX+logM(2);  % polyval with newX
    [~, indexforclosestFrameRt] = min(abs(newX - log10(frameRate))); 
    logdiffusionConstant = newYline(indexforclosestFrameRt);
    diffusionConstant_loggraph = 10^logdiffusionConstant;
    
%     % To visualize the 60th frame MSD value
%     plot(newX,newYline,'g-', 'LineWidth', 2); hold on; 
%     plot(newX(indexforclosestFrameRt), newYline(indexforclosestFrameRt), 'bo', 'MarkerFaceColor', [0 0 0])

    D_log = diffusionConstant_loggraph/4; % MSD = 2dD(dT)^alpha, the final step is to divide by 2d. JD added 9/27/23 
    alpha = logM(1); % Alpha is the slope
end

% Formats the linear graph -- makes the axis in milliseconds, not frames
function [numFrames] = formatLineargraph(durationSec, numFrames)
    xt = xticks; %this is in frames
    xt_ms = (xt).*durationSec./numFrames.*1000; %converted from frames to seconds to milliseconds
    decimals_round = 0;
    [labels_x] = getLabels(xt_ms, decimals_round);
    xticklabels(labels_x{1});
end

% Formats the log10 graph -- converts the log10 axes to 10^x (or 10^y) and also makes x-axis milliseconds, not frames
function [posD, posAlpha] = formatLoggraph(durationSec, numFrames)
    xt = xticks; % This is in frames
    xt_ms = (10.^(xt)).*durationSec./numFrames.*1000; %converted from frames to seconds to milliseconds
    decimals_round = 0;
    [labels_x] = getLabels(xt_ms, decimals_round);
    xticklabels(labels_x{1});
 
    yt = yticks;
    yt_ms = (10.^(yt));
    decimals_round = 1; % Sometimes make this 2 (e.g. very caged/low D/low alpha conditions)
    [labels_y] = getLabels(yt_ms, decimals_round);
    yticklabels(labels_y{1})
    
    posD = [(xt(1)+xt(2))/2 yt(end-1)]; % Position on graph to place the D text
    posAlpha = [(xt(1)+xt(2))/2 (yt(end-1)+yt(end-2))/2]; % Position on graph to place the alpha text
end

% Reformats tick-labels so that I can set them as new xticklabel or yticklabel 
function [labels_xy] = getLabels(xyt_ms, decimal_round)
    eachtick = round(xyt_ms, decimal_round);
    
    labelsString = [];
    for i = 1:length(eachtick)
        labelsString = [labelsString {num2str(eachtick(i))}];
    end
    
    labels_xy{1} = labelsString; % Saving array of new labels for the axis
end