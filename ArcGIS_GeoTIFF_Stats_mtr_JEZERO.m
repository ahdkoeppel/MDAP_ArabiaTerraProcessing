%% Code for Reading in, plotting and comparing GeoTiffs from ArcGIS
% Ari Koeppel - 2019

%Notes: if doing this again, it would be useful to do robust
%filtering/denoising first so PCA doesn't pick out bad params.

clear
close all;


%% Options
Dir = uigetdir('X:\akoeppel\Jezero'); %gets directory
Cluster_option = 1; % 0 = Assign to base clusters only, 1 = assign sophisticated mixture values, 2= assign simple mixture values, XXX = NaN ambiguous pixels
Outliers = 1; %1 = do not remove any outliers. 2 = remove outliers
DOI_name = ["TI"]; %Data of interest
%OtherTI = ["Crater"];
clusters = 10;
P = 1; % percent to subsample
incl_Alb = 0; %1 = Include albedo in PCA, 0 = do not
decorr_slope = 0; %1 = remove variables that are correlated with R1330, 0 = don't

%% Import CRISM Band Parameter Names
opts = detectImportOptions('Y:\Arabia_Terra_MDAP\Statistics\CRISM_L_Params.txt');
CRISM_L_Params = table2array(readtable('Y:\Arabia_Terra_MDAP\Statistics\CRISM_L_Params.txt',opts));
opts = detectImportOptions('Y:\Arabia_Terra_MDAP\Statistics\CRISM_MTR_Params.txt');
CRISM_MTR_Params = table2array(readtable('Y:\Arabia_Terra_MDAP\Statistics\CRISM_MTR_Params.txt',opts));
opts = detectImportOptions('Y:\Arabia_Terra_MDAP\Statistics\CRISM_ratioed_Params.txt');
CRISM_ratioed_Params = table2array(readtable('Y:\Arabia_Terra_MDAP\Statistics\CRISM_ratioed_Params.txt',opts));
opts = detectImportOptions('Y:\Arabia_Terra_MDAP\Statistics\CRISM_AK2014_Params.txt');
CRISM_AK2014_Params = table2array(readtable('Y:\Arabia_Terra_MDAP\Statistics\CRISM_AK2014_Params.txt',opts));
opts = detectImportOptions('Y:\Arabia_Terra_MDAP\Statistics\CRISM_AK2014_Params_New.txt');
CRISM_AK2014_Params_New = table2array(readtable('Y:\Arabia_Terra_MDAP\Statistics\CRISM_AK2014_Params_New.txt',opts));
%% Import Data
Files = dir(fullfile(Dir,'*.tif')); %gets all files
% fullFileName = fullfile(Dir, Files.name);
token = 0;
firstFileName = fullfile(Dir, Files(1).name);
for k = 1:length(Files)
  fullFileName = fullfile(Dir, Files(k).name);
  names = string(Files(k).name(1:end-4));
  fprintf(1, 'Now reading %s\n', fullFileName);
  %opts = detectImportOptions(fullFileName);
  %[RAW,R] = readgeoraster(fullFileName);
  [RAW,R] = geotiffread(fullFileName);
  RAW = double(RAW);
  if k == 1
      if size(RAW,3)>1
          token = 1;
          NumBands = size(RAW,3);
          for i=1:size(RAW,3)
            Band = RAW(:,:,i);
            Bandname = sprintf('set%0.0f_b%0.0f_',token,i);
            if i == 1
            Data = array2table(Band(:),'VariableNames',strcat(Bandname,names));
            else
            Data = horzcat(Data,array2table(Band(:),'VariableNames',strcat(Bandname,names)));
            end
          end
      else
          Data = array2table(RAW(:),'VariableNames',names);
      end
  else
      if size(RAW,3)>1 && token == 0
          token = 1;
          NumBands = size(RAW,3);
          for i=1:size(RAW,3)
            Band = RAW(:,:,i);
            Bandname = sprintf('set%0.0f_b%0.0f_',token,i);
            Data = horzcat(Data,array2table(Band(:),'VariableNames',strcat(Bandname,names)));
          end
      elseif size(RAW,3)>1 && token > 0
          token = 1+token;
          NumBands = size(RAW,3);
          for i=1:size(RAW,3)
            Band = RAW(:,:,i);
            Bandname = sprintf('set%0.0f_b%0.0f_',token,i);
            Data = horzcat(Data,array2table(Band(:),'VariableNames',strcat(Bandname,names)));
          end
      elseif token == 0    
          Data = horzcat(Data,array2table(RAW(:),'VariableNames',names));
      else
          Data = horzcat(Data,array2table(RAW(:),'VariableNames',names));
      end
   end
end

% Rename Vars
for j = 1:width(Data)
    if NumBands == 46 %isempty(strfind(Data.Properties.VariableNames{j},'TRR3')) ~= 1 || isempty(strfind(Data.Properties.VariableNames{j},'trr3')) ~= 1
        num = 1;
        while num <= 46
            strf = sprintf('set%0.0f_b%0.0f_',token,num);
            if isempty(strfind(Data.Properties.VariableNames{j},strf)) ~= 1
                Data.Properties.VariableNames{j} = CRISM_L_Params{num};
                num = num+1;
            end
            if num <= 46 && j < width(Data)
                j = j + 1;
            end
        end
        b=46;
        Albedo = find(contains(Data.Properties.VariableNames,'R1330'));
        Albedo_Name = ['R1330'];
    elseif NumBands >= 60 %isempty(strfind(Data.Properties.VariableNames{j},'mtr')) ~= 1 || isempty(strfind(Data.Properties.VariableNames{j},'MTR')) ~= 1
        num = 1;
        while num <= 60
            strf = sprintf('set%0.0f_b%0.0f_',token,num);
            if isempty(strfind(Data.Properties.VariableNames{j},strf)) ~= 1
                Data.Properties.VariableNames{j} = CRISM_MTR_Params{num};
                num = num+1;
            end
            if num <= 60 && j < width(Data)
                j = j+1;
            end
        end
        b = 60;
        Albedo = find(contains(Data.Properties.VariableNames,'R1330'));
        Albedo_Name = ['R1330'];
    elseif NumBands == 15 %isempty(strfind(Data.Properties.VariableNames{j},'ratioed_Sample')) ~= 1 || isempty(strfind(Data.Properties.VariableNames{j},'lulu')) ~= 1
        for num = 1:15
            strf = sprintf('set%0.0f_b%0.0f_',token,num);
            if isempty(strfind(Data.Properties.VariableNames{j},strf)) ~= 1
                Data.Properties.VariableNames{j} = CRISM_ratioed_Params{num};
            end
            if num < 15 && j < width(Data)
                j = j+1;
            end
        end
        b = 15;
         Albedo = find(contains(Data.Properties.VariableNames,'BD148'));
        Albedo_Name = ['BD148'];
    elseif NumBands == 24 %isempty(find(contains(Data.Properties.VariableNames,'cube_2014'))) ~= 1 && isempty(find(contains(Data.Properties.VariableNames,'Albedo'))) == 1
        for t = 1:token
            for num = 1:24
                strf = sprintf('set%0.0f_b%0.0f_',t,num);
                if isempty(strfind(Data.Properties.VariableNames{j},strf)) ~= 1
                    ParamString = sprintf('_stamp_%0.0f',t);
                    Data.Properties.VariableNames{j} = strcat(CRISM_AK2014_Params_New{num},ParamString);
                end
                if num < 24 && j < width(Data)
                    j = j+1;
                end
            end
        end
        b = 24;
        Albedo = find(contains(Data.Properties.VariableNames,'R1330'));
        Albedo_Name = ['R1330'];
    elseif NumBands == 23 %isempty(find(contains(Data.Properties.VariableNames,'cube_2014'))) ~= 1 && isempty(find(contains(Data.Properties.VariableNames,'Albedo'))) == 0
        for num = 1:23
            strf = sprintf('set%0.0f_b%0.0f_',token,num);
            if isempty(strfind(Data.Properties.VariableNames{j},strf)) ~= 1
                Data.Properties.VariableNames{j} = CRISM_AK2014_Params{num};
            end
            if num < 23 && j < width(Data)
                j = j+1;
            end
        end
        b = 23;
        Albedo = find(contains(Data.Properties.VariableNames,'Albedo'));
        if isempty(Albedo) == 1
            Albedo = find(contains(Data.Properties.VariableNames,'ISLOPE1'));
        end
        Albedo_Name = Data.Properties.VariableNames{Albedo(1)};
    end
    if exist('Albedo','var') == 0
        continue
    elseif isempty(Albedo) == 1
        Albedo = find(contains(Data.Properties.VariableNames,'Albedo'));
        Albedo_Name = Data.Properties.VariableNames{Albedo(1)};
    end
    if j==NumBands||j==NumBands+1||j==NumBands+2
        break
    end
end

%Remove bad variable and NaN no-value pixels
k=width(Data);
i = 1;
while i < k
    Data{Data{:,i}<-10^8,i} = NaN;
    Data{Data{:,i}==65535,i} = NaN;
    Data{Data{:,i}==-65535,i} = NaN;
    Data{Data{:,i}>10^6,i} = NaN;
    if isnumeric(Data{1000,i}) ~= 1 || sum(Data{1:100,i})==0 || all(isnan(Data{:,i}))
        rmvar = string(Data.Properties.VariableNames{i});
        Data(:,rmvar) = []; %remove that variable's column
        i = i-1;
        k = k-1;
    end
    i = i+1;
end
CTX = find(contains(Data.Properties.VariableNames,'CTX'));
Data{Data{:,CTX}==255,CTX} = NaN;
Data{Data{:,end}<-10^8,end} = NaN;
DOI = find(contains(Data.Properties.VariableNames,DOI_name),1);
Data{Data{:,DOI}<=0,DOI} = NaN;
if isnumeric(Data{1000,end}) ~= 1 || sum(Data{1:100,end})==0 || all(isnan(Data{:,end}))
        rmvar = string(Data.Properties.VariableNames{end});
        Data(:,rmvar) = []; %remove that variable's column
end
Data{Data{:,end}<-10^8,end} = NaN;
Data{Data{:,end}==65535,end} = NaN;
Data{Data{:,end}==-65535,end} = NaN;
Data{Data{:,end}>10^6,end} = NaN;


% Find CRISM start and end bands
if b == 46
    startband = find(contains(Data.Properties.VariableNames,CRISM_L_Params{1}));
    endband = find(contains(Data.Properties.VariableNames,CRISM_L_Params{46}));
    while isempty(endband)
        b = b-1;
        endband = find(contains(Data.Properties.VariableNames,CRISM_L_Params{b}));
    end
elseif b >= 60
    startband = find(contains(Data.Properties.VariableNames,CRISM_MTR_Params{1}));
    endband = find(contains(Data.Properties.VariableNames,CRISM_MTR_Params{60}));
    b=60;
    while isempty(endband)
        b = b-1;
        endband = find(contains(Data.Properties.VariableNames,CRISM_MTR_Params{b}));
    end
elseif b == 15
    startband = find(contains(Data.Properties.VariableNames,CRISM_ratioed_Params{1}));
    endband = find(contains(Data.Properties.VariableNames,CRISM_ratioed_Params{15}));
    while isempty(endband)
        b = b-1;
        endband = find(contains(Data.Properties.VariableNames,CRISM_ratioed_Params{b}));
    end
elseif b == 24
    startband = find(contains(Data.Properties.VariableNames,CRISM_AK2014_Params_New{1}));
    endband = find(contains(Data.Properties.VariableNames,CRISM_AK2014_Params_New{24}));
    while isempty(endband)
        b = b-1;
        endband = find(contains(Data.Properties.VariableNames,CRISM_AK2014_Params_New{b}));
    end
elseif b == 23
    startband = find(contains(Data.Properties.VariableNames,CRISM_AK2014_Params{1}));
    endband = find(contains(Data.Properties.VariableNames,CRISM_AK2014_Params{23}));
    while isempty(endband)
        b = b-1;
        endband = find(contains(Data.Properties.VariableNames,CRISM_AK2014_Params{b}));
    end
else
    error('ERROR: The bandnames didn''t print out properly.');
end
Data_orig = Data;

i = startband;
while i <= endband
    if (Data.Properties.VariableNames{i}(1) == 'R' || Data.Properties.VariableNames{i}(1) == 'I' ...
            || isequal(Data.Properties.VariableNames{i}(1:3),'VAR') || isequal(Data.Properties.VariableNames{i}(1:2),'SH') || isequal(Data.Properties.VariableNames{i},'BD1435')) && ~isequal(Data.Properties.VariableNames{i},Albedo_Name)
        Data(:,i) = [];
        endband = endband-1;
    else
        i = i+1;
    end
end
%% Data of Interest
names = Data.Properties.VariableNames;
CRISM = find(contains(Data.Properties.VariableNames,'SINDEX2'),1);
DOI = find(contains(Data.Properties.VariableNames,DOI_name),1);
%DOI2 = find(contains(Data.Properties.VariableNames,OtherTI),1);
Albedo = find(contains(Data.Properties.VariableNames,Albedo_Name),1);


Params1 = ["BD2500_2","D2300","OLINDEX3","LCPINDEX2","BD1900_2",Albedo_Name];
Params2 = ["BD2500_2","D2300","LCPINDEX2","BD1900_2",Albedo_Name,"OLINDEX3"];
Params3 = ["LCPINDEX2","OLINDEX3","D2300", "BD2500_2"];%["R1330","OLINDEX3","BD530_2", "BD860_2"];

if isempty(DOI)==1 || isempty(CRISM) ==1 || isempty(Albedo) == 1
    error('ERROR: Key variable is missing or wrongly labelled.');
end

%% Correlation Calculations Pre Clustering
ScaledData = zeros(height(Data),width(Data));
NewData = zeros(height(Data),width(Data));
for i = 1:width(Data)
    ScaledData(:,i) = rescale(Data{:,i});
    NewData(:,i) = Data{:,i};
end
[CC,~,RL,RU] = corrcoef(ScaledData(:,:),'Rows','complete','Alpha',0.05);

%% PCA
% USEFUL FOR IDing most important variables for classification/clustering
% Do this to remove CRISM outliers, CHECK to see if changes things too much
% Can remove outliers from a specifically noisy band by setting i = to band
% number in data -- NEEDS to BE fine TUNED:
% [B,TF] = rmoutliers(ScaledData(:,5),'percentiles',[0.1 99.5]); %'mean','median','grubbs','gesd'
% ScaledData(TF,5)=NaN;
% ScaledData(:,5) = rescale(ScaledData(:,5));
clear coeff_sorted
clear CSInd
clear impt_vars
clear ScaledData_forPCA
clear ScaledData
for i = 1:width(Data)
    ScaledData(:,i) = rescale(Data{:,i});
    NewData(:,i) = Data{:,i};
end
for i = startband:endband
    if Outliers > 1
        TF = isoutlier(ScaledData(:,i),'grubbs'); %'mean','median','grubbs','gesd'
        % mean = 3 stdev (99.7%) from mean
        ScaledData(TF,i)=NaN;
        ScaledData(:,i) = rescale(ScaledData(:,i));
    end
end

% if incl_Alb >= 1
%     if max(Albedo) <= max(endband) && min(Albedo) >= min(startband)
%         if DOI == 1 
%             ScaledData_forPCA=ScaledData(:,[DOI,min(startband):max(endband)]);
%             PCA_names = names{[DOI,min(startband):max(endband)]};
%         elseif startband == 1
%             ScaledData_forPCA=ScaledData(:,[min(startband):max(endband),DOI]);
%             PCA_names = names{[min(startband):max(endband),DOI]};
%         else
%             ScaledData_forPCA=ScaledData;
%             PCA_names = names;
%         end
%     else
%         ScaledData_forPCA=ScaledData;
%         PCA_names = names;
% %         ScaledData_forPCA=ScaledData(:,[min(startband):max(endband),DOI,Albedo]); %For PCA and clustering, use only physically important variables
% %         PCA_names = names{[min(startband):max(endband),DOI,Albedo]};
%     end
% else
    ScaledData_forPCA=ScaledData;
    PCA_names = names;
% end

if decorr_slope == 1
    Slope_Vars = find(abs(CC(Albedo,:))>0.5); %picks out which variables are strongly correlated to albedo/aspect/topography
    %Slope_Vars = find(or(abs(CC(Albedo,:))>0.45,abs(CC(12,:))>0.45));
    Slope_Vars(Slope_Vars==DOI) = []; %removes TI from this list if it's found to correlate
    Slope_Vars(Slope_Vars==CRISM) = [];
    PCA_names(:,Slope_Vars)=[];
    firstPC = 1; %use first PC below in picking out important vars
    ScaledData_forPCA(:,Slope_Vars)=[];
% elseif sum(explained(2:end)) > 55
%     firstPC = 2; %skip first PC below in picking out important vars b/c it's likely dictated by topography
else
    firstPC = 1;
end
    
% Only use ScaledData_forPCA in calculating important variables and clustering-- But need to plot full variableset

[coeff,~,~,~,explained,mu] = pca(ScaledData_forPCA);

%CV = cov(ScaledData(:,:),'omitrows');
i = firstPC;
while i == 1 || sum(explained(firstPC:i)) <= 55 % use just enough PCs to explain most of the variance
    [coeff_sorted(:,i), CSInd(:,i)] = sort(coeff(:,i)); % sort important PCs
    if i == firstPC
        impt_vars = CSInd(abs(coeff_sorted(:,i))>.25,firstPC); % select highly-weighted variables for each PC (0.3 is a guess)
    else
        impt_vars = cat(1,impt_vars,CSInd(abs(coeff_sorted(:,i))>.25,i)); %Get the indexes of the most important variables for each addtl PC
    end
    i = i+1;
end
impt_vars = unique(impt_vars);

%can visualize PCs with PC1 = coeff(:,1)'*ScaledData_forPCA'; CC_cluster_index = PC1';

%% Gaussian Mixture Model
% options = statset('Display','final','MaxIter',500); 
% GMModel = fitgmdist(ScaledData_forPCA(:,impt_vars),clusters,'Options',options,'CovarianceType','diagonal');
clusters = 10;

AIC = zeros(1,clusters);
GMModels = cell(1,clusters);
options = statset('MaxIter',500);
for k = 1:clusters
    GMModels{k} = fitgmdist(ScaledData_forPCA(:,impt_vars),k,'Options',options,'CovarianceType','diagonal');
    AIC(k)= GMModels{k}.AIC;
end
[minAIC,numComponents] = min(AIC);
GMModel = GMModels{numComponents};
clusters = numComponents;
colors = colormap(hsv(clusters));

GMMidx = cluster(GMModel,ScaledData_forPCA(:,impt_vars));
GMMpost = posterior(GMModel,ScaledData_forPCA(:,impt_vars)); 

figure(1)
Plots = 4;
for i = 1:Plots %width(Data)-1
    subplot(Plots,1,i); %subplot(Plots/2,2,ij);
    ii = find(contains(Data.Properties.VariableNames,Params3{i}));
    for j = 1:clusters
        scatter(ScaledData(GMMidx==j,ii),NewData(GMMidx==j,DOI),5,'o','MarkerEdgeColor','none','MarkerFaceColor',colors(j,:),'MarkerFaceAlpha',1)
        hold on
    end
    hold off
    ylabel(names(DOI),'Interpreter','none')
    xlabel(names(ii),'Interpreter','none')
    legend('1','2','3','4','5','6','Location','best');
end
sgtitle('GMM Hard Clustering')
% clustermarks = ["+","x","o","."];
% for i = 1:Plots %width(Data)-1
%     subplot(Plots,1,i); %subplot(Plots/2,2,ij);
%     ii = find(contains(Data.Properties.VariableNames,Params3{i}));
%     for j = 1:clusters
%         scatter(ScaledData(GMMidx==j,ii),NewData(GMMidx==j,DOI),5,post(GMMidx==j,2),clustermarks{j});
%         hold on
%     end
%     hold off
%     ylabel(names(DOI),'Interpreter','none')
%     xlabel(names(ii),'Interpreter','none')
%     legend('1','2','3','4','Location','best');
% end
% sgtitle('GMM Hard Clustering')
% clrmap = jet(80);
% colormap(clrmap(9:72,:))
% ylabel(colorbar,'Component 3 Posterior Probability')

CC_cluster_index = GMMidx;
%% Reassign Clusters based on likelihood
CC_cluster_index = GMMidx;
%   MC = zeros(length(CC_cluster_index),1); %Max clusters

% Define mixtures of clusters with unique ID #s corresponding to each combo
if Cluster_option == 1
    for i= 1:length(CC_cluster_index)
      [M,In] = max(GMMpost(i,:));
      if isnan(CC_cluster_index(i))
          continue
      elseif M > 0.8 %M = post(i,In)
          CC_cluster_index(i) = In;
    %           MC(i) = 1;
          continue
      elseif sum(GMMpost(i,:)>M/4)==2 %single max is 0.8
          [Result,LocResult] = ismember(find(GMMpost(i,:)>M/4),nchoosek(1:clusters,2),'rows');
          CC_cluster_index(i)=clusters+LocResult;
    %           MC(i) = 1;
          continue
      elseif sum(GMMpost(i,:)>M/4)==3 %single max is 0.66
          %[Result,LocResult] = ismember(find(post(i,:)>M/4),nchoosek(1:clusters,3),'rows');
          CC_cluster_index(i)=clusters+length(nchoosek(1:clusters,2))+1; %could just do +1 here, don't care as much about which
    %       CC_cluster_index(i)=clusters+length(nchoosek(1:clusters,2))+LocResult;%could just do +1 here, don't care as much about which
          continue
      elseif sum(GMMpost(i,:)>M/4)==4 %single max is 0.57
          %[Result,LocResult] = ismember(find(post(i,:)>M/4),nchoosek(1:clusters,4),'rows');
          CC_cluster_index(i)=clusters+length(nchoosek(1:clusters,2))+2;
    % CC_cluster_index(i)=clusters+length(nchoosek(1:clusters,2))+length(nchoosek(1:clusters,3))+LocResult; %could just do +2 here, don't care as much about which
          continue
      else
          CC_cluster_index(i) = clusters+length(nchoosek(1:clusters,2))+2; %NaN;
      end
    end
elseif Cluster_option > 1
    for i= 1:length(CC_cluster_index)
      [M,In] = max(GMMpost(i,:));
      if isnan(CC_cluster_index(i))
          continue
      elseif M > 0.9 %M = post(i,In)
          CC_cluster_index(i) = In;
          continue
      elseif sum(GMMpost(i,:)>0.3)==2 %single max is 0.8
          [Result,LocResult] = ismember(find(GMMpost(i,:)>0.3),nchoosek(1:clusters,2),'rows');
          CC_cluster_index(i)=clusters+LocResult;
          continue
      elseif sum(GMMpost(i,:)>0.2)==3 %single max is 0.66
%           [Result,LocResult] = ismember(find(post(i,:)>0.2),nchoosek(1:clusters,3),'rows');
          CC_cluster_index(i)=clusters+length(nchoosek(1:clusters,2))+1; %could just do +1 here, don't care as much about which
          continue
      elseif sum(GMMpost(i,:)>0.2)==4 %single max is 0.57
%           [Result,LocResult] = ismember(find(post(i,:)>0.2),nchoosek(1:clusters,4),'rows');
          CC_cluster_index(i)=clusters+length(nchoosek(1:clusters,2))+2; %could just do +2 here, don't care as much about which
          continue
      else %if Cluster_option > 2
          CC_cluster_index(i) = clusters+length(nchoosek(1:clusters,2))+2; %NaN;
      end
    end
end

%% Correlation Calculations Within Clusters
Clust_CC=zeros(size(CC,1),size(CC,2),clusters);
Clust_RL=Clust_CC;
Clust_RU=Clust_CC;
% R_squared=Clust_CC;
% CC_cluster_index(CC_cluster_index==4)= 2;
colors = colormap(hsv(clusters));
for i = 1:clusters
    [Clust_CC(:,:,i),~,Clust_RL(:,:,i),Clust_RU(:,:,i)] = corrcoef(ScaledData(CC_cluster_index==i,:),'Rows','complete','Alpha',0.05);
%     mdl = fitlm(ScaledData(CC_cluster_index==i,:),'Rows','complete','Alpha',0.05);
%     R_squared(:,:,i) = mdl.Rsquared.Ordinary;
end

%randomization when subsampling
[m,n] = size(ScaledData); 
% rng('default'); %Use same random seed for reproducibility
% knidx = randperm(m);
% PlottingX = ScaledData(knidx(1:round(P*m)),:);
% PlottingY = NewData(knidx(1:round(P*m)),:);
% CC_cluster_index = CC_cluster_index(knidx(1:round(P*m)));

%Do not Randomize when clustering!
PlottingX = ScaledData;
PlottingY = NewData;

Plots = 4;
figure(2)
for i = 1:Plots %width(Data)-1
    subplot(Plots,1,i); %subplot(Plots/2,2,ij);
    ii = find(contains(Data.Properties.VariableNames,Params3{i}));
    hold on
    for j = 1:clusters
        scatter(PlottingY(CC_cluster_index==j,ii),PlottingY(CC_cluster_index==j,DOI),10,'o','filled','MarkerFaceColor',colors(j,:),'MarkerFaceAlpha',1)
%         RegMod = fitlm(PlottingX(CC_cluster_index==j,ii),PlottingX(CC_cluster_index==j,DOI));
        str{:,j} = sprintf('CC%0.0f: %0.4f +/- %0.4f (@95%%)',j,Clust_CC(ii,DOI,j), Clust_RU(ii,DOI,j)-Clust_CC(ii,DOI,j));
    end
    hold off
    ylabel(names(DOI),'Interpreter','none')
    xlabel(names(ii),'Interpreter','none')
    if clusters == 10
        legend(str{:,1},str{:,2},str{:,3},str{:,4},str{:,5},str{:,6},str{:,7},str{:,8},str{:,9},str{:,10},'Location','best','Interpreter','none');
    elseif clusters == 7
        legend(str{:,1},str{:,2},str{:,3},str{:,4},str{:,5},str{:,6},str{:,7},'Location','best','Interpreter','none');
    elseif clusters == 5
        legend(str{:,1},str{:,2},str{:,3},str{:,4},str{:,5},'Location','best','Interpreter','none');
    elseif clusters == 4
        legend(str{:,1},str{:,2},str{:,3},str{:,4},'Location','best','Interpreter','none');
    elseif clusters == 3
        legend(str{:,1},str{:,2},str{:,3},'Location','best','Interpreter','none'); 
    end
end
sgtitle('Cluster Correlations')

% figure(3)
% for i = 1:Plots %width(Data)-1
%     subplot(Plots,1,i); %subplot(Plots/2,2,ij);
%     ii = find(contains(Data.Properties.VariableNames,Params3{i}));
%     hold on
%     for j = 1:clusters
%         scatter(PlottingY(CC_cluster_index==j,DOI),PlottingY(CC_cluster_index==j,ii),10,'o','filled','MarkerFaceColor',colors(j,:),'MarkerFaceAlpha',1)
%         RegMod = fitlm(PlottingX(CC_cluster_index==j,DOI),PlottingX(CC_cluster_index==j,ii));
%         str{:,j} = sprintf('CC%0.0f: %0.4f, R^2 = %0.4f',j,Clust_CC(DOI,ii,j),RegMod.Rsquared.Ordinary);
%     end
%     hold off
%     xlabel(names(DOI),'Interpreter','none')
%     ylabel(names(ii),'Interpreter','none')
%     if clusters == 10
%         legend(str{:,1},str{:,2},str{:,3},str{:,4},str{:,5},str{:,6},str{:,7},str{:,8},str{:,9},str{:,10},'Location','best','Interpreter','none');
%     elseif clusters == 5
%         legend(str{:,1},str{:,2},str{:,3},str{:,4},str{:,5},'Location','best','Interpreter','none');
%     elseif clusters == 4
%         legend(str{:,1},str{:,2},str{:,3},str{:,4},'Location','best','Interpreter','none');
%     elseif clusters == 3
%         legend(str{:,1},str{:,2},str{:,3},'Location','best','Interpreter','none'); 
%     end
% end
% sgtitle('Cluster Correlations (TI x-axis)')
%% Write new TIFF with Clusters
  [A,R] = geotiffread(firstFileName);
  info = geotiffinfo(firstFileName); %iminfo
  Cluster_File = reshape(CC_cluster_index,size(Band));  %Something may be wrong with the dimensions here
  newfilename = [Files(1).folder,'\Clusters_Mixed',fullFileName(end-3:end)];
  geotiffwrite(newfilename,Cluster_File,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)
  
  out1 = [Files(1).folder,'\GMMidx_mixed.mat'];
  out2 = [Files(1).folder,'\PCAnames.mat'];
  out3 = [Files(1).folder,'\PCA_variables.mat'];
  out4 = [Files(1).folder,'\GMMpost.mat'];
  save(out1,'GMMidx')
  save(out2,'PCA_names')
  save(out3,'impt_vars')
  save(out4,'GMMpost')
  %visualize PCs using PC1 = coeff(:,1)'*ScaledData_forPCA'; CC_cluster_index = PC1';

%% Extract Characteristic TI values
for c = 1:clusters
    Mean(c) = mean(Data{CC_cluster_index==c,DOI}); %Means
    Median(c) = median(Data{CC_cluster_index==c,DOI});
    SD(c) = std(Data{CC_cluster_index==c,DOI}); %standard deviations
    N(c) = length(Data{CC_cluster_index==c,DOI}); %N
end
varNames = ["Mean","Median","SD", "N"];
TI_Stats = table(Mean',Median',SD',N','VariableNames',varNames);
out5 = [Files(1).folder,'\TI_Stats.mat'];
save(out5,'TI_Stats')
% sprintf('Cluster 1 TI: %0.2f +- %0.2f',C1_TI,C1_SD)
% sprintf('Cluster 2 TI: %0.2f +- %0.2f',C2_TI,C2_SD)
% sprintf('Cluster 3 TI: %0.2f +- %0.2f',C3_TI,C3_SD)
% sprintf('Cluster 4 TI: %0.2f +- %0.2f',C4_TI,C4_SD)
% sprintf('Cluster 5 TI: %0.2f +- %0.2f',C5_TI,C5_SD)

% TI_SDs = [C1_TI C2_TI C3_TI C4_TI C5_TI; C1_SD C2_SD C3_SD C4_SD C5_SD]'

% ind = or(CC_cluster_index==1,CC_cluster_index==3);
% C13_TI = mean(Data{ind,DOI});
% C13_SD = std(Data{ind,DOI});
% 
% C12_TI = mean(Data{or(CC_cluster_index==1,CC_cluster_index==2),DOI});
% C12_SD = std(Data{or(CC_cluster_index==1,CC_cluster_index==2),DOI});

%% True axis Cluster Figures
figure(5)
colors = colormap(hsv(clusters));
Plots = 6; %# of plots to show
% The ii = find... fxn has problems when multiple vlaues (ie 2 TI) are
% returned & scatter can't plot
for i = 1:Plots %width(Data)-1
    subplot(Plots/2,2,i); %(width(Data)-1)/2,2,i);
    ii = find(contains(Data.Properties.VariableNames,Params2{i}));
%     clsts1 = gscatter(Data{:,ii},Data{:,DOI},idx,'bgm','.'); *'o'
%     clsts1(1,:).MarkerFaceColor = 'b';clsts1(2,:).MarkerFaceColor = 'g';clsts1(3,:).MarkerFaceColor = 'm';
    for j = 1:clusters
        scatter(Data{CC_cluster_index==j,ii},Data{CC_cluster_index==j,DOI},5,'o','MarkerEdgeColor','none','MarkerFaceColor',colors(j,:),'MarkerFaceAlpha',0.1)
        hold on
    end
    hold off
    ylabel(names(DOI),'Interpreter','none')
    xlabel(names(ii),'Interpreter','none')
end

%% Initial Scaled Density Plots
[m,n] = size(NewData);
%Randomize when subsampling!

%Use same random seed for reproducibility
% knidx = randperm(m);
% P = 1; 
% rng('default'); %Use same random seed for reproducibility
% knidx = randperm(m);
% ScaledData = ScaledData(knidx(1:round(P*m)),:);
% NewData = NewData(knidx(1:round(P*m)),:);

%Do not Randomize when clustering!


figure(1)
Plots = 1; %# of plots to show
for i = 1:Plots %width(Data)-1
    %subplot(Plots/2,2,i); %(width(Data)-1)/2,2,i);
    ii = find(contains(Data.Properties.VariableNames,Params1{i}),1);
    %index = zeros(length(Data{:,ii}),1);
    NewData(isnan(ScaledData(:,DOI)),ii) = NaN;
    [h1, counts] = dscatter(NewData(~isnan(ScaledData(:,ii)),ii),NewData(~isnan(ScaledData(:,ii)),DOI));
    cb = colorbar();
    cb.Ticks = [0,1];
    cb.TickLabels = [1 max(counts)];
    cb.Label.String = 'Counts';
    str1 = sprintf('CC: %0.4f +%0.4f -%0.4f (@95%%)', CC(ii,DOI), RU(ii,DOI)-CC(ii,DOI),CC(ii,DOI)-RL(ii,DOI));
% 	text(mean(ScaledData(:,ii),'omitnan'),mean(ScaledData(:,DOI),'omitnan'),{str},'FontSize',8,'FontName','Ariel');
    title({str1},'FontSize',8,'FontName','Ariel');
    ylabel(names(DOI),'Interpreter','none')
    xlabel(names(ii),'Interpreter','none')
end
sgtitle('Unscaled Density Plots 1')
figure(2)
for ij = 1:Plots
    %subplot(Plots/2,2,ij); %subplot(Plots/2,2,ij);
    ik = find(contains(Data.Properties.VariableNames,Params2{ij}));
    ikm = max(ik);
    %index2 = zeros(length(ScaledData(:,ik)),1);
    NewData(isnan(ScaledData(:,CRISM)),ikm) = NaN;
    [h2, counts2] = dscatter(NewData(~isnan(ScaledData(:,ikm)),ikm),NewData(~isnan(ScaledData(:,ikm)),CRISM));
    cb = colorbar();
    cb.Ticks = [0,1];
    cb.TickLabels = [1, max(counts2)];
    cb.Label.String = 'Counts';
    str2 = sprintf('CC: %0.4f +%0.4f -%0.4f (@95%%)', CC(ikm,CRISM), RU(ikm,CRISM)-CC(ikm,CRISM),CC(ikm,CRISM)-RL(ikm,CRISM));
%     text(mean(ScaledData(:,ik),'omitnan'),mean(ScaledData(:,DOI),'omitnan'),{str2},'FontSize',8,'FontName','Ariel');
    title({str2},'FontSize',8,'FontName','Ariel');
    ylabel(names(CRISM(1)),'Interpreter','none')
    xlabel(names(ikm),'Interpreter','none')
end
sgtitle('Unscaled Density Plots 2')
%Use these plots to estimate the number of clusters

%% Functions

function [hAxes, counts] = dscatter(X,Y, varargin)
% DSCATTER creates a scatter plot coloured by density.
%
%   DSCATTER(X,Y) creates a scatterplot of X and Y at the locations
%   specified by the vectors X and Y (which must be the same size), colored
%   by the density of the points.
%
%   DSCATTER(...,'Marker',M) allows you to set the marker for the
%   scatter plot. Default is 'o', circle.
%
%   DSCATTER(...,'Msize',MS) allows you to set the marker size for the
%   scatter plot. Default is 10.
%
%   DSCATTER(...,'Filled',false) sets the markers in the scatter plot to be
%   outline. The default is to use filled markers.
%
%   DSCATTER(...,'Plottype',TYPE) allows you to create other ways of
%   plotting the scatter data. Options are "surf','mesh' and 'contour'.
%   These create surf, mesh and contour plots colored by density of the
%   scatter data.
%
%   DSCATTER(...,'Bins',[NX,NY]) allows you to set the number of bins used
%   for the 2D histogram used to estimate the density. The default is to
%   use the number of unique values in X and Y up to a maximum of 200.
%
%   DSCATTER(...,'Smoothing',LAMBDA) allows you to set the smoothing factor
%   used by the density estimator. The default value is 20 which roughly
%   means that the smoothing is over 20 bins around a given point.
%
%   DSCATTER(...,'Logy',true) uses a log scale for the yaxis.
%
%   Examples:
%
%       [data, params] = fcsread('SampleFACS');
%       dscatter(data(:,1),10.^(data(:,2)/256),'log',1)
%       % Add contours
%       hold on
%       dscatter(data(:,1),10.^(data(:,2)/256),'log',1,'plottype','contour')
%       hold off
%       xlabel(params(1).LongName); ylabel(params(2).LongName);
%       
%   See also FCSREAD, SCATTER.
% Copyright 2003-2004 The MathWorks, Inc.
% $Revision:  $   $Date:  $
% Reference:
% Paul H. C. Eilers and Jelle J. Goeman
% Enhancing scatterplots with smoothed densities
% Bioinformatics, Mar 2004; 20: 623 - 628.
Lambda = [];
Nbins = [];
Plottype = 'scatter';
ContourFlag = false;
Msize = 5;
Marker = 'o';
Logy = false;
Filled = true;
if nargin > 2
    if rem(nargin,2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'Smoothing','Bins','Plottype','Logy','Marker','Msize','Filled'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs); %#ok
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % smoothing factor
                    if isnumeric(pval)
                        Lambda = pval;
                    else
                        error('Bioinfo:InvalidScoringMatrix','Invalid smoothing parameter.');
                    end
                case 2
                    if isscalar(pval)
                        Nbins = [ pval pval];
                    else
                        Nbins = pval;
                    end
                case 3
                    Plottype = pval;
                case 4
                    Logy = pval;
                    Y = log10(Y);
                case 5
                    ContourFlag = pval;
                case 6
                    Marker = pval;
                case 7
                    Msize = pval;
                case 8
                    Filled = pval;
            end
        end
    end
end
minx = min(X,[],1);
maxx = max(X,[],1);
miny = min(Y,[],1);
maxy = max(Y,[],1);
if isempty(Nbins)
    Nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200) ];
end
if isempty(Lambda)
    Lambda = 20;
end
edges1 = linspace(minx, maxx, Nbins(1)+1);
ctrs1 = edges1(1:end-1) + .5*diff(edges1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(miny, maxy, Nbins(2)+1);
ctrs2 = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];
[n,p] = size(X);
bin = zeros(n,2);
% Reverse the columns to put the first column of X along the horizontal
% axis, the second along the vertical.
[dum,bin(:,2)] = histc(X,edges1);
cts = max(dum);
[dum,bin(:,1)] = histc(Y,edges2);
cts(2) = max(dum);
counts = max(cts);
H = accumarray(bin,1,Nbins([2 1])) ./ n;
G = smooth1D(H,Nbins(2)/Lambda);
F = smooth1D(G',Nbins(1)/Lambda)';
% = filter2D(H,lambda);
if Logy
    ctrs2 = 10.^ctrs2;
    Y = 10.^Y;
end
okTypes = {'surf','mesh','contour','image','scatter'};
k = strmatch(lower(Plottype), okTypes); %#ok
if isempty(k)
    error('dscatter:UnknownPlotType',...
        'Unknown plot type: %s.',Plottype);
elseif length(k)>1
    error('dscatter:AmbiguousPlotType',...
        'Ambiguous plot type: %s.',Plottype);
else
    switch(k)
        case 1 %'surf'
            h = surf(ctrs1,ctrs2,F,'edgealpha',0);
        case 2 % 'mesh'
            h = mesh(ctrs1,ctrs2,F);
        case 3 %'contour'
            [dummy, h] =contour(ctrs1,ctrs2,F);
        case 4 %'image'
            nc = 256;
            F = F./max(F(:));
            colormap(repmat(linspace(1,0,nc)',1,3));
            h =image(ctrs1,ctrs2,floor(nc.*F) + 1);
        case 5 %'scatter'
            F = F./max(F(:));
            ind = sub2ind(size(F),bin(:,1),bin(:,2));
            col = F(ind);
            if Filled
                h = scatter(X,Y,Msize,col,Marker,'filled');
            else
                h = scatter(X,Y,Msize,col,Marker);
            end
    end
end
if Logy
    set(gca,'yscale','log');
end
if nargout > 0
    hAxes = get(h,'parent');
end
end

function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
end

function [outliers, DBclusts] = DBOPT(eps,minpts,Data)
    DBinx = dbscan(Data,eps,minpts,'Distance','squaredeuclidean');
    outliers = -sum(DBinx==-1);
    DBclusts = -max(DBinx);
end
