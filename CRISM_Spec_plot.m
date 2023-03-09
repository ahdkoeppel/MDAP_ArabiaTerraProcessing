%% Code for Reading in and plotting CRISM spectra
% Ari Koeppel - 2019
clear
close all;

%% Import Data
FILE = char('Y:\Arabia_Terra_MDAP\Spectra\W_of_Yelapa\FRT0000B710_trr3_Clays.txt');
opts = detectImportOptions(FILE);
Tbl = readtable(FILE,opts);
count = 0;
for i = 1:size(Tbl,2)
    if isnumeric(Tbl{1,i}) && count == 0
        l = i;
        count = 1;
    elseif isnumeric(Tbl{1,i}) && count == 1
        s = i;
    end
end        
lambda = Tbl{:,l}; 
if lambda(1)>10
    lambda = lambda./1000;
end
spectra = Tbl{:,s};
for i = 1:length(spectra)
    if spectra(i) > 2 || spectra(i) < 0.1
        spectra(i) = NaN;
    end
    if i > 1 && abs(spectra(i)-spectra(i-1)) > 0.1
        spectra(i) = NaN;
    end
    
end
    
%% Smoothing
smooth_spec = smoothdata(spectra,'sgolay');
%'(r)lowess' — Linear regression over each window of A. This method can be computationally expensive, but results in fewer discontinuities.
%'sgolay' — Savitzky-Golay filter, which smooths according to a quadratic polynomial that is fitted over each window of A. This method can be more effective than other methods when the data varies rapidly.
%'(r)loess' — Robust quadratic regression over each window of A. This method is a more computationally expensive version of the method 'loess', but it is more robust to outliers.
%'gaussian' — Gaussian-weighted moving average over each window of A.
%% Plot
plot(lambda,spectra);
%plot(lambda,zeros(length(lambda),1));
xlabel('Wavelength (\mum)','FontName','Arial');
ylabel('(Ratioed) Relative Reflectance','FontName','Arial')
% ind1 = strfind(FILE,'\frt');
% if isempty(ind1) 
%     ind1 = strfind(FILE,'\frs');
% end
% if isempty(ind1) 
%     ind1 = strfind(FILE,'\hrl');
% end
% if isempty(ind1) 
%     ind1 = strfind(FILE,'\hrs');
% end
% ind2 = strfind(FILE,'.txt');
% title(FILE(ind1(end)+1:ind2(end)-1),'Interpreter','none');
title(FILE(1:end-4),'Interpreter','none');
ax = gca; ax.XLim = [1.0145 2.6];

%CO2
% line('XData',[1.44 1.44],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.58 1.58],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.99 1.99],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.29 2.29],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.35 2.35],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');

%H2O
% line('XData',[1.03 1.03],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.26 1.26],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.48 1.48],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.02 2.02],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .5 .8],'LineWidth',.5,'LineStyle','--');

% %HCP
% line('XData',[1.015 1.015],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .8 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.16 2.16],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .8 .8],'LineWidth',.5,'LineStyle','--');
% 
% %Ol
% line('XData',[1.12 1.12],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .9 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.3 1.3],'YData',get(gca, 'ylim'),...
%     'Color',[.8 .9 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.05 1.05],'YData',get(gca, 'ylim'),...
%     'Color',[.9 .8 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.25 1.25],'YData',get(gca, 'ylim'),...
%     'Color',[.9 .8 .8],'LineWidth',.5,'LineStyle','--');

% %PH Sulfate
% line('XData',[1.43 1.43],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.94 1.94],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.43 2.43],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .8],'LineWidth',.5,'LineStyle','--');
% 
% %MH Sulfate
% line('XData',[1.63 1.63],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .5 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.97 1.97],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .5 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.14 2.14],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .5 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.40 2.40],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .5 .1],'LineWidth',.5,'LineStyle','--');
% 
% %Mg Smectite
% line('XData',[1.41 1.41],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.92 1.92],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.31 2.31],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.39 2.39],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .8 .4],'LineWidth',.5,'LineStyle','--');
% 
% %Fe/Mg/Ca Carbonate
% line('XData',[1.91 1.91],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .4 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.34 2.34],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .4 .8],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.53 2.53],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .4 .8],'LineWidth',.5,'LineStyle','--');
% 
% %Bassanite
% line('XData',[1.43 1.43],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .2 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.77 1.77],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .2 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.92 1.92],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .2 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.26 2.26],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .2 .4],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.49 2.49],'YData',get(gca, 'ylim'),...
%     'Color',[.4 .2 .4],'LineWidth',.5,'LineStyle','--');
% 
% %Gypsum
% line('XData',[1.44 1.44],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.49 1.49],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.54 1.54],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.75 1.75],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.95 1.95],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.21 2.21],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.26 2.26],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.48 2.48],'YData',get(gca, 'ylim'),...
%     'Color',[.1 .2 .1],'LineWidth',.5,'LineStyle','--');
% 
% %OH,Fe,Sulfate
% line('XData',[1.48 1.48],'YData',get(gca, 'ylim'),...
%     'Color',[.7 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.82 1.82],'YData',get(gca, 'ylim'),...
%     'Color',[.7 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.99 1.99],'YData',get(gca, 'ylim'),...
%     'Color',[.7 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.19 2.19],'YData',get(gca, 'ylim'),...
%     'Color',[.7 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.23 2.23],'YData',get(gca, 'ylim'),...
%     'Color',[.7 .2 .1],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.36 2.36],'YData',get(gca, 'ylim'),...
%     'Color',[.7 .2 .1],'LineWidth',.5,'LineStyle','--');
% 
% %Jarosite
% line('XData',[1.47 1.47],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.85 1.85],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.27 2.27],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');

% %Muscovite
% line('XData',[1.41 1.41],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.93 1.93],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.21 2.21],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.35 2.35],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.45 2.45],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');

%Al-Smectite
% line('XData',[1.41 1.41],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[1.91 1.91],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');
% line('XData',[2.20 2.20],'YData',get(gca, 'ylim'),...
%     'Color',[.2 .2 .7],'LineWidth',.5,'LineStyle','--');

