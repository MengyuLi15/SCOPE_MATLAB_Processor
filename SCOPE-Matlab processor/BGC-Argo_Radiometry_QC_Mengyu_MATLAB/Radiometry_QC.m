%% 14/01/2025 SCOPE
%
% Take 1902605 as example
% MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
% Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
% method via: Organelli et al., 2016 JTECH, DOI: 10.1175/JTECH-D-15-0193.1
%
%
%    QC=1 good data
%    QC=2 probably good
%    QC=3 probably bad 
%    QC=4 bad
%
% !!! this matlab code need several extra function in the mean file dir
%
% including:
%
%       SolarAzEl.m
%       BGCArgo_single_prof_QC.m
%       RadiometryProfilePointsNumTest.m
%       RadiometryFirstPolynomialFitTest.m
%       RadiometrySecondPolynomialFitTest.m
%       RTQCRadiometry.m
%       CalculationKdFromEd.m
%       CalculationZeuZpdFromPARProfile.m
%       pop_kw.m

clc;clear;close all;
files = dir('B*.nc');
Argo_REFERENCE_DATE_TIME = datetime('1950-01-01 00:00:00'); % juld date reference
lat_all = nan(1,1); lon_all = nan(1,1); UTC_DATE_all = nan(1,1);
Kd380_all = nan(1,1); Kd412_all = nan(1,1); Kd490_all = nan(1,1);
Zeu_all =  nan(1,1); flag_Zeu_Chla = nan(1,1);

for i = 133:length(files)

filename = files(i).name;% have the file name, it is okay to loop later
% ncdisp(filename) % to check the ncinfo
warning('off');

PRES = ncread(filename,'PRES');
LATITUDE = double(ncread(filename,'LATITUDE'));% location
LONGITUDE = double(ncread(filename,'LONGITUDE'));% location
TIME = ncread(filename,'JULD');% time
UTC_DATE = TIME+juliandate(Argo_REFERENCE_DATE_TIME);
UTC_DATE = datetime(UTC_DATE,'convertfrom','juliandate');% utc time
LATITUDE = LATITUDE(1,1);LONGITUDE = LONGITUDE(1,1);UTC_DATE = UTC_DATE(1,1);

IRR_380 = ncread(filename,'DOWN_IRRADIANCE380');% W/m^2/nm
IRR_380_QC = ncread(filename,'DOWN_IRRADIANCE380_QC');
[IRR_380,press_IRR_380] = BGCArgo_single_prof_QC(PRES,IRR_380,IRR_380_QC);
IRR_380 = IRR_380.*100; % convert to uW cm-2 nm-1
% semilogy(press_IRR_380,IRR_380)
IRR_412 = ncread(filename,'DOWN_IRRADIANCE412');% W/m^2/nm
IRR_412_QC = ncread(filename,'DOWN_IRRADIANCE412_QC');
[IRR_412,press_IRR_412] = BGCArgo_single_prof_QC(PRES,IRR_412,IRR_412_QC);
IRR_412 = IRR_412.*100; % convert to uW cm-2 nm-1

IRR_490 = ncread(filename,'DOWN_IRRADIANCE490');% W/m^2/nm
IRR_490_QC = ncread(filename,'DOWN_IRRADIANCE490_QC');
[IRR_490,press_IRR_490] = BGCArgo_single_prof_QC(PRES,IRR_490,IRR_490_QC);
IRR_490 = IRR_490.*100; % convert to uW cm-2 nm-1

PAR = ncread(filename,'DOWNWELLING_PAR');% microMoleQuanta/m^2/sec
DOWNWELLING_PAR_QC = ncread(filename,'DOWNWELLING_PAR_QC');
[PAR,press_PAR] = BGCArgo_single_prof_QC(PRES,PAR,DOWNWELLING_PAR_QC);
flag_Zeu_Chla_t = nan;
% in the function BGCArgo_single_prof_QC.m, I have already removed the
% 99999 and nan in the profile

% initialization QC flag variable
flag_IRR_380 = nan(length(IRR_380),1);
flag_IRR_412 = nan(length(IRR_412),1);
flag_IRR_490 = nan(length(IRR_490),1);
flag_PAR = nan(length(PAR),1);

flag_IRR_380 = RTQCRadiometry(IRR_380,press_IRR_380,flag_IRR_380,...
    UTC_DATE,LATITUDE,LONGITUDE);
flag_IRR_412 = RTQCRadiometry(IRR_412,press_IRR_412,flag_IRR_412,...
    UTC_DATE,LATITUDE,LONGITUDE);
flag_IRR_490 = RTQCRadiometry(IRR_490,press_IRR_490,flag_IRR_490,...
    UTC_DATE,LATITUDE,LONGITUDE);
flag_PAR = RTQCRadiometry(PAR,press_PAR,flag_PAR,...
    UTC_DATE,LATITUDE,LONGITUDE);

%% calculate kd and zeu
% first check the zeu and zpd
[Zeu,Zpd] = CalculationZeuZpdFromPARProfile(PAR,press_PAR,flag_PAR);
% if there is no PAR profile, we can use the Chla profile and their relationship 
% to calculate the Zeu
%
if isnan(Zeu)
    try
        Chla = ncread(filename,'CHLA_ADJUSTED');
        Chla_QC = ncread(filename,'CHLA_ADJUSTED_QC');
        [Chla,press_Chla] = BGCArgo_single_prof_QC(PRES,Chla,Chla_QC);
        [~,b] = min(abs(press_Chla-5));
        chla_5 = Chla(b);
        Zeu = 34.0.*chla_5.^(0-0.39); % Lee et al., (2007), JGR Oceans,https://doi.org/10.1029/2006JC003802
        Zpd = Zeu./4.6;
        % disp(['Zeu is from Chla! Filename is: ',filename])
        flag_Zeu_Chla_t = 1;
        if isempty(Zeu)
            Zeu = nan; Zpd = nan;
        end
    catch ME
        if (strcmp(ME.identifier,'MATLAB:imagesci:netcdf:unknownLocation'))
            Chla = [];
            Chla_QC = [];
            Zeu = nan;
            flag_Zeu_Chla_t = nan;
            % disp(['Zeu is nan! Filename is: ',filename])
        end
    end
end

% try
% Chla = ncread(filename,'CHLA_ADJUSTED');
% Chla_QC = ncread(filename,'CHLA_ADJUSTED_QC');
% catch ME
%     if (strcmp(ME.identifier,'MATLAB:imagesci:netcdf:unknownLocation'))
%         Chla = [];
%         Chla_QC = [];
%         Zeu = nan;
%         flag_Zeu_Chla_t = nan;
%     else
%         [Chla,press_Chla] = BGCArgo_single_prof_QC(PRES,Chla,Chla_QC);
%         if isnan(Zeu)
%         [~,b] = min(abs(press_Chla-5));
%         chla_5 = Chla(b);
%         Zeu = 34.0.*chla_5.^(0-0.39); % Lee et al., (2007), JGR Oceans,https://doi.org/10.1029/2006JC003802
%         Zpd = Zeu./4.6;
%         disp(['Zeu is from Chla! Filename is: ',filename])
%         flag_Zeu_Chla_t = 1;
%         else
%             flag_Zeu_Chla_t = 0;
%         end
%     end
% end


Zeu_all = [Zeu_all;Zeu];
%% calculate Kd from Ed
% [Kd_zpd,SE_zpd] = CalculationKdFromEd(Ed,press_Ed,flag_Ed,Zpd,lambda)
if isempty(Zpd) | isnan(Zpd)
    Kd380_zpd = nan; SE380_zpd = nan;
    Kd412_zpd = nan; SE412_zpd = nan;
    Kd490_zpd = nan; SE490_zpd = nan;
    flag_Zeu_Chla_t = nan;
else
[Kd380_zpd,SE380_zpd] = CalculationKdFromEd(IRR_380,...
    press_IRR_380,flag_IRR_380,Zpd,380);
[Kd412_zpd,SE412_zpd] = CalculationKdFromEd(IRR_412,...
    press_IRR_412,flag_IRR_412,Zpd,412);
[Kd490_zpd,SE490_zpd] = CalculationKdFromEd(IRR_490,...
    press_IRR_490,flag_IRR_490,Zpd,490);
if isnan(flag_Zeu_Chla_t)
    flag_Zeu_Chla_t = 0;
end
end
%%
Kd380_all = [Kd380_all;Kd380_zpd];
Kd412_all = [Kd412_all;Kd412_zpd];
Kd490_all = [Kd490_all;Kd490_zpd];
lat_all = [lat_all;LATITUDE];
lon_all = [lon_all;LONGITUDE];
UTC_DATE_all = [UTC_DATE_all;datenum(UTC_DATE)];
flag_Zeu_Chla = [flag_Zeu_Chla;flag_Zeu_Chla_t];
%% Save results into each mat file
% filename = [filename(1:14),'_Radiometry.mat'];
% save(filename,'flag*','IRR*','PAR','LATITUDE','LONGITUDE','UTC_DATE','press*','Z*','Kd*');
% 
% % Depth = (1:1:Zpd)';
% % ProfileTime = repmat(UTC_DATE,length(Depth),1);
% % Lon = repmat(LONGITUDE,length(Depth),1);
% % Lat = repmat(LATITUDE,length(Depth),1);
% Output_tablehead = table(UTC_DATE,LONGITUDE,LATITUDE,Zeu,Zpd);
% 
% Output_table = table(Kd380_zpd,SE380_zpd,Kd412_zpd,SE412_zpd,Kd490_zpd,SE490_zpd);
% 
% filename = [filename(1:14),'_Radiometry.xlsx'];
% 
% writetable(Output_tablehead,filename,'Sheet',1,'Range','A1')
% writetable(Output_table,filename,'Sheet',1,'Range','A4')

%% Visualization
% subplot(1,4,1);% 380
% Radiometry_final = IRR_380;
% Radiometry_final_press = press_IRR_380;
% Radiometry_final(flag_IRR_380==3) = [];
% Radiometry_final_press(flag_IRR_380==3) = [];
% 
% semilogx(IRR_380,press_IRR_380,'DisplayName','Measured','MarkerFaceColor',[0.89 0.80 0.80],...
%     'MarkerSize',3,...
%     'Marker','o',...
%     'LineWidth',1,...
%     'Color',[0.80 0.80 0.80]);hold on
% semilogx(Radiometry_final,Radiometry_final_press,'DisplayName','QC','MarkerFaceColor',...
%     [0.851 0.325 0.098],'MarkerSize',2,...
%     'Marker','o',...
%     'LineWidth',1);
% 
% axis('ij');
% xlabel('E_d(380) [µW cm^{-2} nm^{-1}]')
% ylabel('Pressure [dbar]')
% legend1 = legend('show');
% set(legend1,'Location','southeast');
% 
% subplot(1,4,2);% 412
% Radiometry_final = IRR_412;
% Radiometry_final_press = press_IRR_412;
% Radiometry_final(flag_IRR_412==3) = [];
% Radiometry_final_press(flag_IRR_412==3) = [];
% 
% semilogx(IRR_412,press_IRR_412,'DisplayName','Measured','MarkerFaceColor',[0.89 0.80 0.80],...
%     'MarkerSize',3,...
%     'Marker','o',...
%     'LineWidth',1,...
%     'Color',[0.80 0.80 0.80]);hold on
% semilogx(Radiometry_final,Radiometry_final_press,'DisplayName','QC','MarkerFaceColor',...
%     [0.851 0.325 0.098],'MarkerSize',2,...
%     'Marker','o',...
%     'LineWidth',1);
% 
% axis('ij');
% xlabel('E_d(412) [µW cm^{-2} nm^{-1}]')
% ylabel('Pressure [dbar]')
% legend1 = legend('show');
% set(legend1,'Location','southeast');
% 
% subplot(1,4,3);% 490
% Radiometry_final = IRR_490;
% Radiometry_final_press = press_IRR_490;
% Radiometry_final(flag_IRR_490==3) = [];
% Radiometry_final_press(flag_IRR_490==3) = [];
% 
% semilogx(IRR_490,press_IRR_490,'DisplayName','Measured','MarkerFaceColor',[0.89 0.80 0.80],...
%     'MarkerSize',3,...
%     'Marker','o',...
%     'LineWidth',1,...
%     'Color',[0.80 0.80 0.80]);hold on
% semilogx(Radiometry_final,Radiometry_final_press,'DisplayName','QC','MarkerFaceColor',...
%     [0.851 0.325 0.098],'MarkerSize',2,...
%     'Marker','o',...
%     'LineWidth',1);
% 
% axis('ij');
% xlabel('E_d(490) [µW cm^{-2} nm^{-1}]')
% ylabel('Pressure [dbar]')
% legend1 = legend('show');
% set(legend1,'Location','southeast');
% 
% subplot(1,4,4);% PAR
% Radiometry_final = PAR;
% Radiometry_final_press = press_PAR;
% Radiometry_final(flag_PAR==3) = [];
% Radiometry_final_press(flag_PAR==3) = [];
% 
% semilogx(PAR,press_PAR,'DisplayName','Measured','MarkerFaceColor',[0.89 0.80 0.80],...
%     'MarkerSize',3,...
%     'Marker','o',...
%     'LineWidth',1,...
%     'Color',[0.80 0.80 0.80]);hold on
% semilogx(Radiometry_final,Radiometry_final_press,'DisplayName','QC','MarkerFaceColor',...
%     [0.851 0.325 0.098],'MarkerSize',2,...
%     'Marker','o',...
%     'LineWidth',1);
% 
% axis('ij');
% xlabel('PAR [µmol quanta m^{-2} s^{-1}]')% microMoleQuanta/m^2/sec
% ylabel('Pressure [dbar]')
% legend1 = legend('show');
% set(legend1,'Location','southeast');
   
end
result = [Kd380_all Kd412_all Kd490_all];
UTC_DATE_all = datetime(UTC_DATE_all,"ConvertFrom","datenum");
result(1,:) = []; lat_all(1,:) = []; lon_all(1,:) = [];
Zeu_all(1,:) = [];UTC_DATE_all(1,:) = [];
flag_Zeu_Chla(1,:) = [];

load splat
sound(y,Fs)
clearvars -except result UTC_DATE_all lat_all lon_all Zeu_all flag_Zeu_Chla
