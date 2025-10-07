function [CphytoProfile,Press_Cphyto] = PhytoplanktonCarbonBellacicco(bbpProfile,Pressure,Zeu,SampleDate,Location)
%
%    Important information:
%       1. bbp should be after quality controlled profile
%       2. wavelength of bbp is only for 700 now (10/03/2025 version)
%       3. bbpProfile and Pressure should be in the same depth, and size;
%       4. The input of location is char, according to Bellacicco et al.,
%           2019, GRL, 10.1029/2019gl084078, it is also necesssary
%       5. If there will be the updates of the bbpk, it can be directly
%           updated in the 'bbpkBellacicco2019.mat' file
%
%
%    PhytoplanktonCarbonBellacicco detail information:
%
%    SCOPE 10/03/2025
%    in MATLAB R2024b
%
%    for SCOPE pj Cphyto 
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Bellacicco et al., 2019, GRL, doi: 10.1029/2019gl084078
%
%    This function is for Cphyto calculation from bbp700 (Quality controlled!) 
%        of the BGC-Argo profiles, according to Bellacicco et al., 2019 GRL
%
%    Inputs：
%
%    bbpProfile：Quality controlled bbp profile, Nx1, m-1;
%    Pressure：Quality controlled depth same size with bbp, Nx1, m; 
%    Zeu：euphotic depth, 1x1, m;
%    SampleDate：sample date of bbp, in MATLAB datetime or datenum are both ok;
%    Location: char, the float location should be in the list of Bellacicco et al., 2019
%
%    Output:
%
%    CphytoProfile: Cphyto profile, Nx1, mg C m-3
%    Press_Cphyto: Depth of the Cphyto, Nx1, m

%        PS: According to the paper:
%       (i) the surface layer: the layer between sea surface and the first optical depth; 
%       (ii) the euphotic layer: the layer between sea surface and euphotic zone;
%       (iii) the bottom layer: the layer between 950 and 1,000 m.


%% test
% bbpProfile = BBP700(:,307);
% Pressure = press_BBP700(:,307);
% Zeu = Zpl.("Z_0.05chl")(307).*1.5;
% SampleDate = UTC_DATE(307);
% Location = "NASTG";

%% Step 0: Find the coefficients of the bbpk upper and lower layers based on the given location
load('bbpkBellacicco2019.mat');
SF = 16455; % from Bellacicco et al., 2019 Sensors
bbp = bbpProfile;
Press_t = Pressure;

if Location == "ASZ_SIZ"
    bbpk_upper = bbpkeuphotic.ASZ_SIZ;
    bbpk_bottom = bbpkbottom.ASZ_SIZ;
elseif Location == "SA"
    bbpk_upper = bbpkeuphotic.SA;
    bbpk_bottom = bbpkbottom.SA;
elseif Location == "STZ"
    bbpk_upper = bbpkeuphotic.STZ;
    bbpk_bottom = bbpkbottom.STZ;
elseif Location == "PFZ"
    bbpk_upper = bbpkeuphotic.PFZ;
    bbpk_bottom = bbpkbottom.PFZ;
elseif Location == "IOMZ"
    bbpk_upper = bbpkeuphotic.IOMZ;
    bbpk_bottom = bbpkbottom.IOMZ;
elseif Location == "NASPG"
    bbpk_upper = bbpkeuphotic.NASPG;
    bbpk_bottom = bbpkbottom.NASPG;
elseif Location == "NS"
    bbpk_upper = bbpkeuphotic.NS;
    bbpk_bottom = bbpkbottom.NS;
elseif Location == "WMS"
    bbpk_upper = bbpkeuphotic.WMS;
    bbpk_bottom = bbpkbottom.WMS;
elseif Location == "EMS"
    bbpk_upper = bbpkeuphotic.EMS;
    bbpk_bottom = bbpkbottom.EMS;
elseif Location == "IEQ"
    bbpk_upper = bbpkeuphotic.IEQ;
    bbpk_bottom = bbpkbottom.IEQ;
elseif Location == "NASTG"
    bbpk_upper = bbpkeuphotic.NASTG;
    bbpk_bottom = bbpkbottom.NASTG;
elseif Location == "NPSTG"
    bbpk_upper = bbpkeuphotic.NPSTG;
    bbpk_bottom = bbpkbottom.NPSTG;
elseif Location == "PSEW"
    bbpk_upper = bbpkeuphotic.PSEW;
    bbpk_bottom = bbpkbottom.PSEW;
elseif Location == "SASTG"
    bbpk_upper = bbpkeuphotic.SASTG;
    bbpk_bottom = bbpkbottom.SASTG;
elseif Location == "SISTG"
    bbpk_upper = bbpkeuphotic.SISTG;
    bbpk_bottom = bbpkbottom.SISTG;
elseif Location == "SPSTG"
    bbpk_upper = bbpkeuphotic.SPSTG;
    bbpk_bottom = bbpkbottom.SPSTG;
else
    warning('Location input is not valid');
    CphytoProfile = nan(length(bbp),1);
    Press_Cphyto = nan(length(bbp),1);
    return
end

%% Step 1: find the month in the bbpk vector
catelog_bbpk = month(SampleDate);
bbpk_upper_this = bbpk_upper(catelog_bbpk);
bbpk_bottom_this = bbpk_bottom(catelog_bbpk);

if isnan(bbpk_upper_this) || isnan(bbpk_bottom_this)
    warning('There is no bbpk in this time');
    CphytoProfile = nan(length(bbp),1);
    Press_Cphyto = nan(length(bbp),1);
    return
end

%% Step 2: Prepare the bbp data for the Cphyto calculation
Press_t(isnan(bbp)) = [];
bbp(isnan(bbp)) = [];

if isempty(bbp)
    warning('There is no data');
    CphytoProfile = nan(length(bbp),1);
    Press_Cphyto = nan(length(bbp),1);
    return
end

[~,index] = rmoutliers(bbp,'movmedian' ,7);
bbp(index) = nan;

Press_t(isnan(bbp)) = [];
bbp(isnan(bbp)) = [];

if isempty(bbp)
    warning('There is no data');
    CphytoProfile = nan(length(bbp),1);
    Press_Cphyto = nan(length(bbp),1);
    return
end

%% Step 3: bbp upper layer and bottom layer

bbp_upper = bbp;
% if Press_t(end) <= Zeu
%     bbp_upper = bbp;
%     bbp_bottom = nan;
% else
%     bbp_upper = bbp(Press_t <= Zeu);
%     bbp_bottom = bbp(Press_t > Zeu);
% end

%% Step 4: Cphyto
Cphyto_upper = (bbp_upper-bbpk_upper_this).*SF;
% Cphyto_bottom = (bbp_bottom-bbpk_bottom_this).*SF;
Cphyto_upper(Cphyto_upper<0) = 0;
% Cphyto_bottom(Cphyto_bottom<0) = 0;

CphytoProfile = Cphyto_upper;
[~,index] = rmoutliers(CphytoProfile,'movmedian',3);
CphytoProfile(index) = nan;
% CphytoProfile = [Cphyto_upper;Cphyto_bottom];
Press_Cphyto = Press_t;
% plot(Press_t, bbp)
% plot(Press_Cphyto, CphytoProfile)
CphytoProfile(Press_t > 500) = 0;

if max(CphytoProfile,[],'all',"omitnan")>100
    CphytoProfile = nan(length(bbp),1);
    Press_Cphyto = nan(length(bbp),1);
    return
end

end
