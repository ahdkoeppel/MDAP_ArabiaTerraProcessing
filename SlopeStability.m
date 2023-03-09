%% SlopeStability.m
%Slope Stability Bishop Analysis
%Ari Koeppel 2021 - modified from 
clear all
close all

% Obs=xlsread('Soil_Cohesion.xlsx');
% NumObs = length(Obs);
NumObs = 1
rho = 1200:100:2000; %kg/m^3
g = 3.711; %m/s^2
W = rho.*g;
phi = 25:5:60; %degrees
c = zeros(length(W),length(phi),NumObs); %Cohesion Matrix
F = 1; %Factor of Safety in equilibrium

NSlices = 10; %number of slices to use

for i = 1:NumObs
%     Xtop = Obs(i,1); %
%     Xbot = Obs(i,2); %
%     Ytop = Obs(i,3); %
%     Ybot = Obs(i,4); %
    Xtop = 301; %
    Xbot = 106; %
    Ytop = 2448; %
    Ybot = 2586; %
    L = abs(Xtop-Xbot);
    H = abs(Ytop-Ybot);
    dl = L/NSlices; %b
    dh = H/NSlices; %h
    for j = 1:length(W)
        for k = 1:length(phi)
            syms X
            alpha = atand(H/L);
            eqn = 1./sum(W(j).*sind(alpha).*ones(length(NSlices))).*sum((X.*dl+W(j)*tand(phi(k)).*secd(alpha)./(1+tand(alpha).*tand(phi(k))./F)).*ones(length(NSlices))) == F;
            c(j,k,i) = solve(eqn,X,'Real',true);
        end
    end
    cmax(i) = max(c(:,:,i),[],'all')/1000; %kPa
    cmin(i) = min(c(:,:,i),[],'all')/1000; %kPa
end
            
