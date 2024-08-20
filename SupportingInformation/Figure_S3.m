% Figure S3
% Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .

clear; close all; clc;

% Data folder
all_folder = '../SupportingInformation/Data_Figure_S3';
cd(all_folder);

subjectID = '10003L';

startDir = strcat( '', cd,'');
fileNames = dir('*.csv');

allMatches = table();

for i = 1: size(fileNames,1)
    
    % Table with experiment results
    if contains(fileNames(i).name, '_S3_RW2paramCtrlExp.csv') == 1
        allMatches = readtable(fileNames(i).name);
    end

end

% Colors for the 7 conditions tested
color(1,:) = [0      0.4470 0.7410];
color(2,:) = [0.8500 0.3250 0.0980];
color(3,:) = [0.9290 0.6940 0.1250];
color(4,:) = [0.4940 0.1840 0.5560];
color(5,:) = [0.4660 0.6740 0.1880];
color(6,:) = [0.3010 0.7450 0.9330];
color(7,:) = [0.6350 0.0780 0.1840];

% Two Background conditions
numBackgrounds = unique(allMatches.Background); % "1" defines Background-present; "2" defines Background-absent

for b = 1: length(numBackgrounds)
    
    % Plotting results: Diffusion constant vs alpha
    figure('color','w', 'Renderer', 'painters', 'Position', [10 10 600 500]); hold on;
    index = 1;
    curBkgd = allMatches(allMatches.Background == numBackgrounds(b),:);
    plot(curBkgd.TestAlpha,curBkgd.TestD, 'o','Color', [0 0 0], 'MarkerFaceColor', [0 0 0],'markersize',12);
    
    numDConditions = unique(curBkgd.TestD);
    for i = 1: length(numDConditions) 
        
        % Plotting two matches each condition
        curCond = curBkgd(curBkgd.TestD == numDConditions(i),:);
        
        for a = 1: size(curCond,1)
            plot(curCond.MatchedAlpha(a),curCond.MatchedD(a), 'o','Color',  color(index,:), 'MarkerFaceColor',  color(index,:),'markersize',7); 
            plot ([curCond.TestAlpha(a) curCond.MatchedAlpha(a)],[curCond.TestD(a) curCond.MatchedD(a)],'Color', color(index,:), 'linewidth',1);
        end
        
        index= index+1;
    end
     
    title(subjectID);
    xlabel('\alpha'); ylabel('Diffusion Constant [arcmin^2/s]');
    xticks(0.6:0.1:1.8); yticks(0:10:90);
    xlim([0.6 1.8]); ylim([0 90]);
    set(gca, 'FontSize', 18, 'fontname','helvetica')
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
    grid on; box on;    
end