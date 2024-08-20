% Figure 5 
% Written by: Josie D'Angelo 2024. Email: josephine_dangelo@berkeley.edu .

clear; close all; clc;

cd('../Data_Figure_5');

files = dir(pwd);
dirFlags = [files.isdir];
folders = files(dirFlags);
folderNames = {folders(3:end).name};
 
for f = 1:  length(folderNames)
    cd(folderNames{f})
    
    csvNames = dir('*.csv');
    
    figure('color','w'); hold on;  
    for t = 1: length(csvNames)
        trace = readtable(csvNames(t).name); 
        plot(trace.X,trace.Y, 'k'); % units are in [arcmin]
    end
    
    plot(0,0,'+','color','r','markersize',10, 'linewidth',4); % Plotting starting position red cross
    plot([0 4],[-10 -10],'k-','linewidth',1.1) % Label for 4 arcmin
    axisLim = 20; 
    xlim([-axisLim axisLim]); ylim([-axisLim axisLim]);
    xticks(-axisLim:2:axisLim); yticks(-axisLim:2:axisLim);
    grid on; axis square; box on;
    title(['10003L' sprintf(' (%1.f traces)', length(csvNames))])
    
    cd ..
end

