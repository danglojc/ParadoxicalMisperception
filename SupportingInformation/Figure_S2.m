% Figure S2
% Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .

clear; close all; clc;
path = '../SupportingInformation/Data_Figure_S2';
cd(path);

files = dir(pwd);
dirFlags = [files.isdir];
subjects = files(dirFlags);
subjectNames = {subjects(3:end).name};

for s = 1: length(subjectNames)
    cd(subjectNames{s});
    subfiles = dir(pwd);
    subdirFlags = [subfiles.isdir];
    subfolders = subfiles(subdirFlags);
    subfolderNames = {subfolders(3:end).name};
    
    for f = 1: length(subfolderNames)
        cd(subfolderNames{f});
        % Plot all traces for each background and Gain condition
        figure('color','w'); hold on;
        csvFiles = dir('*.csv');
        for t = 1: length(csvFiles)
            trace = readtable(csvFiles(t).name); 
            plot(trace.X,trace.Y, 'k'); 
        end
        
        plot(0,0,'+','color','r','markersize',10, 'linewidth',4); % Plotting starting position red cross
        plot([0 8],[-25 -25],'k-','linewidth',1.1) % Label for 8 arcmin
        axisLim = 30;
        xlim([-axisLim axisLim]); ylim([-axisLim axisLim])
        xticks(-axisLim:5:axisLim); yticks(-axisLim:5:axisLim);
        grid on; axis square; box on;
        title([strrep(subfolderNames{f},'_', ' ') sprintf(' (%1.f traces)', length(csvFiles))])
        
        cd .. % exiting Gain/Bkgd folder
    end
    cd .. % exiting subject's folder
end