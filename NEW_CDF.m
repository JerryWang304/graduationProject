%d1 = dlmread('C:\\b_T0.txt',' ');
%d2 = dlmread('C:\\b_T1.txt',' ');

d1 = dlmread('F:\\ECMP.txt',' ');
d2 = dlmread('F:\\wait_and_hop_time.txt',' ');

filename = 'B_CDF';

linestyle = ['--','-'];
color = ['r','b'];
legendpos = 'Best';
marker = ['s','o'];
lw = 4.0;
ms = 10;
fs = 16;
% legendkey = {'Per Flow','T=1','T=2','T=4','T=6', 'T=12'};
% legendkey = {'N=12','N=6','N=4','N=3','N=2', 'N=1'};
% legendkey = {'Per Flow','T=1','T=5','T=10'};
legendkey = {'ECMP','Wait and Hop'};


[h1, stat1] = cdfplot(d1);
% set(h1, 'LineStyle','-','color','r','LineWidth',lw, 'Marker','o','MarkerSize',4);
set(h1, 'LineStyle','-','color','r','LineWidth',lw);
hold on;
[h2, stat2] = cdfplot(d2);
set(h2, 'LineStyle','-','color','b','LineWidth',lw);

hold off;
title('');
xlabel('Total Completion Time','FontSize', fs, 'FontName', 'Arial');
ylabel('');
legend(legendkey,'Location', 'SouthEast');
set(gcf,'position',[100 100 640 320]);
set(gca, 'FontSize', fs, 'FontName', 'Arial','YGrid','on');
%set(gca, 'xscale', 'log');
set(gcf,'PaperPositionMode','auto');
print('-r0','-depsc', strcat(filename, '.eps'));
% ps2pdf -dEPSCrop result_GP_per.eps
clear;