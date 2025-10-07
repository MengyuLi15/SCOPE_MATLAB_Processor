function [POCProfile,QC_POC] = ParticulateOrganicCarbonKoestner(ChlaProfile,bbpProfile,Pressure,lambda)
%
%    Important information:
%       1. bbp and Chla should be after quality controlled profiles
%       2. wavelength of bbp is necesscary, only 700 and 550
%       3. Chla and bbp should be in the same depth, and size;
%       4. QC of the POC including:
%               < 0 test
%               global range test
%               sigma inf test
%               Chla outliner test (movemedian 7)
%               bbp outliner test (movemedian 7)
%
%    ParticulateOrganicCarbonKoestner detail information:
%
%    SCOPE 17/2/2025
%    in MATLAB R2024b
%
%    for SCOPE pj POC 
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Koestner et al., 2022 FMARS, DOI: 10.3389/fmars.2022.941950
%
%    This function is for POC calculation from bbp700, bbp550 and Chla (Quality controlled!) 
%        of the BGC-Argo profiles, according to Koestner et al., 2022 FMARS
%
%    Inputs：
%
%    ChlaProfile：Quality controlled Chla profile, Nx1, mg m-3;
%    bbpProfile：Quality controlled bbp700 or bbp550 profile, Nx1, m-1; 
%    Pressure：Quality controlled depth, Nx1, m;
%    lambda：bbp wavelength, 1x1, 550 or 700 nm;
%
%    Output:
%
%    POCProfile: POC profile, Nx1, mg m-3
%    QC_POC: QC flag of POC profile
%            QC good           flag = 1
%            QC probebaly good flag = 2 (not used)
%            QC probebaly bad  flag = 3
%            QC bad            flag = 4
%
%
%% test
% ChlaProfile = Chlat;
% bbpProfile = BBP700t;
% Pressure = press_BBP700t;
% lambda = 700;
%% Step 0: Different fitting coefficients are selected according to different bands
QC_POC = nan(length(ChlaProfile),1);
if lambda == 700
    k1 = 89.423; k2 = 0.1881; k3 = 0.7591; k4 = 0.1934;
    e1 = 1.636; e2 = 21.2; POCmin = 33.4;
elseif lambda == 550
    k1 = 206.16; k2 = 0.3615; k3 = 0.6623; k4 = 0.1504;
    e1 = 2.013; e2 = 34.9; POCmin = 34.4;
else
    POCProfile = nan(length(ChlaProfile),1);
    QC_POC(isnan(QC_POC)) = 4;
    warning('POC cannot be calculated without the corresponding coefficient in this wavelength');
    return
end

%% Step 1: Preparation
Chla_temp = ChlaProfile;
bbp_temp = bbpProfile;
press_temp = Pressure;
index_num = (1:1:length(ChlaProfile))';
QC_POC(Chla_temp<0 | bbp_temp<0) = 4;% smaller than 0 flag = 4
bbp_temp(bbp_temp<0) = nan;
Chla_temp(Chla_temp<0) = nan;
QC_POC(Chla_temp>12.15 | bbp_temp > 0.03) = 4; % globel range test, Organelli et al., (2017), ESSD
bbp_temp(bbp_temp > 0.03) = nan;
Chla_temp(Chla_temp>12.15) = nan;
QC_POC(bbp_temp==0) = 3; % sigma = inf
QC_POC(press_temp > 800 & bbp_temp > 0.0004) = 3;

index = isnan(Chla_temp) | isnan(bbp_temp);
QC_POC(index) = 4;
Chla_temp(index) = []; bbp_temp(index) = []; press_temp(index) = [];index_num(index) = [];

if isempty(Chla_temp) | isempty(bbp_temp) | isempty(press_temp)
    POCProfile = nan(length(ChlaProfile),1);
    QC_POC(isnan(QC_POC)) = 4;
    warning('There is no data for POC calculation');
    return   
end
[~,TF] = rmoutliers(Chla_temp,'movmedian',7); % Chla outliners
QC_POC(index_num(TF)) = 3;
[~,TF] = rmoutliers(bbp_temp,'movmedian',7); % bbp outliners
    QC_POC(index_num(TF)) = 3;
%% Step 2: sigma ς calcualtion

sigmaraio = ChlaProfile./bbpProfile; % mg m-2

%% Step 3: POC* calcualtion
%     POC_star = k1.*(bbp700.^k2).*(slope.^k3).*(slope.^(k4.*log10(bbp700)));

POCstar = k1.*(bbpProfile.^k2).*(sigmaraio.^k3).*(sigmaraio.^(k4.*log10(bbpProfile)));

%% Step 4: POC calculation

POCProfile = nan(length(ChlaProfile),1);

for j = 1:length(POCProfile)
    if POCstar(j) < POCmin
        POCProfile(j) = e1.*POCstar(j)-e2;
    elseif POCstar(j) >= POCmin
        POCProfile(j) = POCstar(j);
    end

end

QC_POC(isnan(QC_POC)) = 1;

% plot(POCProfile,Pressure)
% POCProfile(QC_POC==3|QC_POC==4) = nan;
% plot(POCProfile,Pressure)
end


