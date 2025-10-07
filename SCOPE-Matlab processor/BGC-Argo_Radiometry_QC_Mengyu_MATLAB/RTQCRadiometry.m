function FlagRadiometry = RTQCRadiometry(Radiometry,PressRadiometry,FlagOriginalRadiometry,...
    UTC_DATE,LATITUDE,LONGITUDE)
%
%    RTQCRadiometry detail information:
%
%    SCOPE 17/1/2025
%    in MATLAB R2024b
%
%    for SCOPE pj PAR AND IRR RT QC 
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Organelli et al., 2016 JTECH, DOI: 10.1175/JTECH-D-15-0193.1
%
%    This function is for the RT Quality Control for irr and par BGC-Argo
%    profiles, according to Organelli et al., 2016 JTECH
%
%    QC=1 good data
%    QC=2 probably good
%    QC=3 probably bad 
%    QC=4 bad
%
%    Notes:
%    1. ! Please leave the Flag = 3 and 4
%    2. !! This matlab function need several extra function in the mean file dir, including:
%
%       SolarAzEl.m
%       BGCArgo_single_prof_QC.m
%       RadiometryProfilePointsNumTest.m
%       RadiometryFirstPolynomialFitTest.m
%       RadiometrySecondPolynomialFitTest.m
%
%    Inputs：
%
%    Radiometry：Radiometry profile, Nx1, irr or par profile;
%    PressRadiometry：Press of Radiometry profile, Nx1
%    FlagOriginalRadiometry：Initialize the flag, Nx1
%    UTC_DATE：Sample time, in UTC time, 1x1, MATLAB datetime;
%    LATITUDE：latitude, 1x1, double
%    LONGITUDE：longitude, 1x1, double
%
%    Output:
%
%    FlagRadiometry: flag after QC, Nx1



FlagRadiometry = FlagOriginalRadiometry;

Radiometry_temp = Radiometry;
PressRadiometry_temp = PressRadiometry;
index = isnan(Radiometry_temp)|isnan(PressRadiometry_temp);
Radiometry_temp(index) = [];
PressRadiometry_temp(index) = [];

if isempty(Radiometry_temp)|isempty(PressRadiometry_temp)|isnan(datenum(UTC_DATE))|isnan(LATITUDE)|isnan(LONGITUDE)
    FlagRadiometry(isnan(FlagRadiometry)) = 4;
    warning('Profile or location is empty! All flag = 4.');
    return;
end

%% Step 0: Sun elevation detection / negative detection
% 1. Sun elevation detection according to location and time 
[~, El] = SolarAzEl(UTC_DATE,LATITUDE,LONGITUDE,0);
if El<2 % if the sun elevation <2, then all is flag=3
    FlagRadiometry(isnan(FlagRadiometry)) = 3;
    warning('Did not pass the sun elevation detection! All flag = 3.')
    return;
end
% 2. negative detection
FlagRadiometry(Radiometry<=0|isnan(Radiometry)) = 4;

Radiometry_temp = Radiometry;
PressRadiometry_temp = PressRadiometry;
index_all = (1:1:length(Radiometry))';

Radiometry_temp(FlagRadiometry==3|FlagRadiometry==4) = [];
PressRadiometry_temp(FlagRadiometry==3|FlagRadiometry==4) = [];
index_all(FlagRadiometry==3|FlagRadiometry==4) = [];

if length(Radiometry_temp)<=5
    warning('There are not enough valid data! All flag = 4.')
    FlagRadiometry(isnan(FlagRadiometry)) = 4;
    return;
end

%% Step 1.1: DARK IDENTIFICATION 
% Check the data and perform the Lilliefors normality test
% Initializes the Lilliefors normality test result vector

% index_temp = (1:1:length(Radiometry_temp))';
LillieforsRadiometry = nan(length(Radiometry_temp) - 4, 1);

for l = 1:length(LillieforsRadiometry)
    try
        % Extracts the current subarray and performs the Lilliefors test
        currentData = Radiometry_temp(l:end);
        [~, pValue] = lillietest(currentData); % Lilliefors 
        LillieforsRadiometry(l) = pValue; % save p
    catch
        % If there is an error, copy the previous value
        if l > 1
            LillieforsRadiometry(l) = LillieforsRadiometry(l - 1);
        else
            LillieforsRadiometry(l) = NaN; % If it is the first value, set it to NaN
        end
    end
end
% put the flag 3 according to the Lilliefors test results
index_flag_temp = find(abs(LillieforsRadiometry)>=0.01); % have a little bit change than the original one
if ~isempty(index_flag_temp)
    flag_temp_firstdepth = index_flag_temp(1,1)+4;
    FlagRadiometry(index_all(flag_temp_firstdepth):end) = 3;
else
end
%% Step 1.2: Profile points detection
% need more than 5 point for the fit after take all flag 3;
% if it is short than 5 points, then make all flag=3
%
[FlagRadiometry, isLessThan5Points] = RadiometryProfilePointsNumTest(Radiometry,PressRadiometry,FlagRadiometry);
if isLessThan5Points
    warning('There are not enough valid data! All flag = 4.')
    FlagRadiometry(isnan(FlagRadiometry)) = 4;
    return;
end

%% Step 2: Identification of clouds, wave focusing, and spikes - first Polynomial fit test
% outliers among residual values produced with respect to a nonlinear fit on radiometric
% profiles after removal of dark measurements.
%
FlagRadiometry = RadiometryFirstPolynomialFitTest(Radiometry,PressRadiometry,FlagRadiometry);
%% Step 3: Identification of clouds, wave focusing, and spikes - second Polynomial fit test
% outliers among residual values produced with respect to a nonlinear fit on radiometric
% profiles after removal of dark measurements.
%
FlagRadiometry = RadiometrySecondPolynomialFitTest(Radiometry,PressRadiometry,FlagRadiometry);

end