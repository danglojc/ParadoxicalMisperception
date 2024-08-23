% Figure 4
% Written by: Josie D'Angelo 2024.

% Violin Plot from:
% Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
% https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847, This code is released under the terms of the BSD 3-clause license.
% Modified 'Violin.m': Josephine D'Angelo added lines 259-282. 

clear; close all; clc;

% Setting path
current_path = pwd;
addpath(current_path);

if ismac
    cd([current_path, '/Data_Figure_4']);
else
    cd([current_path, '\Data_Figure_4']);
end

individGains= readtable('All_Eye_individTrials.csv');

% Defining parameters
subjects      = unique(individGains.Subject);
numBackground = unique(individGains.Background); % "1" defines Background-present; "2" defines Background-absent
numGains      = unique(individGains.Gain);

meansTable = table();
for s = 1:length(subjects)
    for b = 1:length(numBackground)
        for g = 1: length(numGains)
            data_of_interest = individGains(individGains.Subject == subjects(s) & individGains.Background == numBackground(b) & individGains.Gain == numGains(g),:);
            newArray = [subjects(s) numBackground(b) numGains(g) mean(data_of_interest.D_Eye) mean(data_of_interest.Alpha_Eye) mean(data_of_interest.D_Perceived)];
            meansTable = [meansTable; array2table(newArray, 'VariableNames', {'Subject', 'Background','Gain','D_Eye','Alpha_Eye', 'D_Perceived'})];
        end
    end
end

color_minus = [0 0 0.6];
color_zero  = [0.9290 0.6940 0.1250];
color_plus  = [0 0.35 0];

%% Figure 4A: Plotting diffusion constant for eye motion 

dataCols = [meansTable.D_Eye(meansTable.Background==1 & meansTable.Gain==-1.5), meansTable.D_Eye(meansTable.Background==1 & meansTable.Gain==1.5), meansTable.D_Eye(meansTable.Background==1 & meansTable.Gain==0)...
            meansTable.D_Eye(meansTable.Background==2 & meansTable.Gain==-1.5), meansTable.D_Eye(meansTable.Background==2 & meansTable.Gain==1.5), meansTable.D_Eye(meansTable.Background==2 & meansTable.Gain==0)];

[~] = plotFigure(dataCols, color_minus, color_plus, color_zero);

yticks(0:5:50); ylim([0 50])
ylabel('Diffusion constant for eye motion, D_{EM} [arcmin^2/s]');

%% Figure 4B: Plotting alpha for eye motion 

dataCols = [meansTable.Alpha_Eye(meansTable.Background==1 & meansTable.Gain==-1.5), meansTable.Alpha_Eye(meansTable.Background==1 & meansTable.Gain==1.5), meansTable.Alpha_Eye(meansTable.Background==1 & meansTable.Gain==0)...
            meansTable.Alpha_Eye(meansTable.Background==2 & meansTable.Gain==-1.5), meansTable.Alpha_Eye(meansTable.Background==2 & meansTable.Gain==1.5), meansTable.Alpha_Eye(meansTable.Background==2 & meansTable.Gain==0)];

[~] = plotFigure(dataCols, color_minus, color_plus, color_zero);

yticks(0:0.1:2);ylim([0.8 2])
ylabel('\alpha for eye motion, \alpha_{EM}');

%% Figure 4C: Plotting speed for eye motion 

speedTable= readtable('All_Speed_Eye.csv');

dataCols = [speedTable.Speed_Eye(speedTable.Background==1 & speedTable.Gain==-1.5), speedTable.Speed_Eye(speedTable.Background==1 & speedTable.Gain==1.5), speedTable.Speed_Eye(speedTable.Background==1 & speedTable.Gain==0)...
            speedTable.Speed_Eye(speedTable.Background==2 & speedTable.Gain==-1.5), speedTable.Speed_Eye(speedTable.Background==2 & speedTable.Gain==1.5), speedTable.Speed_Eye(speedTable.Background==2 & speedTable.Gain==0)];

[~] = plotFigure(dataCols, color_minus, color_plus, color_zero);

yticks(0:5:50); ylim([0 50])
ylabel('Speed for eye motion, S_{EM} [arcmin/s]');

cd ..; % Returning to ParadoxicalMisperception folder
%% Function to implement 'Violin.m'
function [dataCols] = plotFigure(dataCols, color_minus, color_plus, color_zero)
    categories = categorical({'-1.5' '+1.5' '0' '-1.5' '+1.5' '0'}); 

    figure('color','w'); hold on;

    % Gain -1.5: Background-present
    vs1 = Violin({dataCols(:,1)},1, "ViolinColor",{color_minus});

    % Gain +1.5: Background-present
    vs2 = Violin({dataCols(:,2)},2, "ViolinColor",{color_plus});

    % Gain 0: Background-present
    vs3 = Violin({dataCols(:,3)},3, "ViolinColor",{color_zero});

    % Gain -1.5: Background-absent
    vs4 = Violin({dataCols(:,4)},4, "ViolinColor",{color_minus});

    % Gain +1.5: Background-absent
    vs5 = Violin({dataCols(:,5)},5, "ViolinColor",{color_plus});

    % Gain 0: Background-absent
    vs6 = Violin({dataCols(:,6)},6, "ViolinColor",{color_zero});

    set(gca,xticklabels=categories)
    xlim([0.5, 6.5]);
    
    set(gca, 'FontSize', 18, 'fontname','helvetica');
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
end
