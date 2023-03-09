%% Script to derive erosion rates from crater size frequency distributions by using estimated production and erosion functions
% Ari Koeppel - 2020
%Adapted from Smith et al., 2008 and Kite et al., 2017

clear all
%close all

[file,path] = uigetfile('X:\akoeppel\Arabia_Terra_MDAP\CraterDensities\DensityMaps\Plotting\*.csv');
raw = readmatrix([path,file]);

%% Read in list of crater diameters in km
Data = raw(:,2)./1000;
%Enter area of measurment in km^2
prompt = 'Enter unit area (km^2): ';
Area = input(prompt);

NumBins = 10;
%[N_bin,Edge_bin] = histcounts(sqrt(2).*log10(Data),NumBins); %Log spaced binning. Could also do sqrt(2) spaced bins (Hartmann et al., 2005)
% for i = 1:NumBins
%     D_bin(i) = (10^(Edge_bin(i))/sqrt(2)+10^Edge_bin(i+1)/sqrt(2))/2;
% end
Edge_bin = [0.0078125;0.011049;0.015625;0.022097;0.03125;0.044194;0.0625;0.088388;0.125;0.17678;0.25;0.35355;0.5;0.70711;1;1.4142;2;2.83]; %sqrt(2) spacing
%Edge_bin = [0.02;0.0251;0.0316;0.0398;0.0501;0.0631;0.0794;0.1;0.126;0.158;0.2;0.251;0.316;0.398;0.501;0.631;0.794;1;1.26;1.58;2;2.51;3.16]; %Michael 2013 10/decade binning
[N_bin] = histcounts(Data,Edge_bin);
D_bin = Edge_bin(1:end-1);
% D_bin = [0.02 0.05 0.1 0.2 0.4 0.8 1 1.5];
% N_bin = [100 90 80 60 40 30 10 2];
for i = 1:length(N_bin)
    N_cum(i) = sum(N_bin(i:end));
end

%% Plot check vs CraterstatsII

figure(1)
scatter(D_bin,N_bin./Area,'s');set(gca,'xscale','log');set(gca,'yscale','log')
title('crater density by size(sqrt(2) bins)')
xlabel('Diameter (km)');ylabel('craters/km^2');
figure(2)
scatter(D_bin,N_cum./Area,'s');set(gca,'xscale','log');set(gca,'yscale','log')
title('cumulative crater density greater than size (sqrt(2) bins)')
xlabel('Diameter (km)');ylabel('craters/km^2');

%% Bin crater sizes based on X system
%X = psudolog, sqrt(2) etc.

%bins = number of craters in each bin
%bindiam = size of crater at center of each bin
%% Production & Erosion Functions Smith 2008
for i = 1:length(N_bin)
    if D_bin(i) >= 0.001 && D_bin(i) < 1.4
        P(i) = 0.0035*(0.13*log(D_bin(i))+0.83)/(D_bin(i)^3.3);
    elseif D_bin(i) >= 1.4 && D_bin(i) < 48.1
        P(i) = 10^(-1.8*log10(D_bin(i))-2.59);
    elseif D_bin(i) >= 48.1
        P(i) = 10^(-2.2*log10(D_bin(i))-1.89);
    end
    if D_bin(i) < 1.2
        phi(i) = 1+0.1225/(0.1225+D_bin(i)^2);
    elseif D_bin(i) >= 1.2
        phi(i) = 0;
    end
    if D_bin(i) < 5.8
        eta(i) = 0.2*D_bin(i);
    elseif D_bin(i) >= 5.8
        eta(i) = (0.42*log(D_bin(i))-0.01);
    end
    %lambda(i) = phi(i)*B/(1000*eta(i));
end


%% Solve for erosion rate
% syms X
% time = 1; %Extract from large craters using crater tools (in Ga)
% for i = 1:length(N_bin)
%     eqns(i) = P(i)/(phi(i)*X/(1000*eta(i)))*(1-exp(-(phi(i)*X/(1000*eta(i)))*time)) == N_bin(i);
%     S(i) = solve(eqns(i),X,'Real',true) %bin-wise estimate of erosion
% end
% solve(eqns,[X],'Real',true)


% syms X time
% for i = 1:length(bins)
%     eqns(i) = P(i)/(phi(i)*X/(1000*eta(i)))*(1-exp(-(phi(i)*X/(1000*eta(i)))*time)) == bins(i);
%     S(i) = solve(eqns(i),[X time],'Real',true);
% end
% solve(eqns,[X time],'Real',true)

%% Kite 2017 
%Step 1: See if site is characterized by steady exhumation (this can be
%done by visual inspection too) -- Maybe skip this step and say assuming
%steady exhumation
D_min = 0.1; % minimum diameter at incompleteness turnoff
alpha = 1+sum(N_bin)*sum(log(D_bin./D_min));
sigma = (alpha-1)/sqrt(sum(N_bin));

%Step 1 Alternate: Plot in CratertoolsII and see if there is a steady fall
%off the isochron

%Step 2: Calculate erosion rate by bin
phi_1 = 0.5;
B = sqrt(2);
for i = 1:length(N_bin)
    %Cumulative Crater flux (get new one from Michael et al 2013) at sqrt(2) intervals
    %starting from 0.00781
    if i == 1
        P(i) = 1.89*10^3;
    elseif i==2
        P(i) = 7.54*10^2;
    elseif i==3
        P(i) = 2.96*10^2;
    elseif i==4
        P(i) = 1.05*10^2;
    elseif i==5
        P(i) = 3.86*10^1;
    elseif i==6
        P(i) = 1.46*10^1;
    elseif i==7
        P(i) = 5.17*10^0;
    elseif i==8
        P(i) = 1.87*10^0;
    elseif i==9
        P(i) = 6.51*10^-1;
    elseif i==10
        P(i) = 2.14*10^-1;
    elseif i==11
        P(i) = 6.66*10^-2;
    elseif i==12
        P(i) = 1.96*10^-2;
    elseif i==13
        P(i) = 5.75*10^-3;
    elseif i==14
        P(i) = 1.73*10^-3;
    elseif i==15
        P(i) = 5.84*10^-4;
    elseif i==16
        P(i) = 2.76*10^-4;
    elseif i ==17
        P(i) = 1.48*10^-4;
    end
%     Michael 2013 correction - doesn't seem right
%     if i ==17
%         k = (P(i)-7.92*10^-5)/(D_bin(i)-2.83);
%     else
%         k = (P(i)-P(i+1))/(D_bin(i)-D_bin(i+1));
%     end
%     P(i) = P(i)*(B^(k/2)-B^(-k/2))/(k*(sqrt(B)-1/sqrt(B)));
    E(i) = 0.2*D_bin(i)*phi_1/(N_cum(i)/P(i)/Area)*1000;
end
figure(3)
loglog(D_bin,E);
xlabel('Diameter (km)');ylabel('Erosion Rate (nm/yr)');
hold on