function [Kd_zpd,SE_zpd] = CalculationKdFromEd(Ed,press_Ed,flag_Ed,Zpd,lambda)
%
%    CalculationKdFromEd detail information:
%
%    SCOPE 20/1/2025
%    in MATLAB R2024b
%
%    for SCOPE pj
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    This function is for calculate Kd within pd
%
%    Inputs：
%    Ed：Ed profile, Nx1, and don t have to QC in advance, uW cm-2 nm-1;
%    press_Ed：Press of Ed profile, Nx1, m
%    flag_Ed：Ed QC flag, Nx1
%    Zpd: euphotic depth Zeu, m
%    lambda: wavelength, nm
%
%    Output:
%    Kd_zpd: kd 0-Zpd, m-1
%    SE_zpd: Standard Error of Kd of 0-Zpd
%
% 
% Zpd = Zeu./4.6;
%% test
% [Kd380_zpd,SE380_zpd] = CalculationKdFromEd(IRR_380,...
%     press_IRR_380,flag_IRR_380,Zpd,380);

% Ed = IRR490zeu; press_Ed = p_IRR490zeu; flag_Ed = ones(length(IRR490zeu),1); Zpd = zeut; lambda = 490;
%% 
Ed_QCed = Ed;
press_Ed_QCed = press_Ed;
Ed_QCed(flag_Ed==3|flag_Ed==4)  = nan;
press_Ed_QCed(isnan(Ed_QCed)) = [];
Ed_QCed(isnan(Ed_QCed)) = [];% after the qc and then for the next calculation

if isempty(Ed_QCed)
    Kd_zpd = nan;
    SE_zpd = nan;
    return
end

% press_bin = (0.5:1:max(press_Ed_QCed))'; % 1 m bin
% Ed_QCed_bin = interp1(press_Ed_QCed,Ed_QCed,press_bin,"pchip");
% Kd_Profile = nan(length(press_bin)-1,1);
% for a = 1:(length(press_bin)-1)
%     Kd_Profile(a) = 1./1.*log(Ed_QCed_bin(a)./Ed_QCed_bin(a+1));
% end
% Press_Kd_Profile = (1:1:press_bin(end)-0.5)';
% 
% % Kd_inZeu
% x = press_bin(press_bin<=Zeu);
% y = log(Ed_QCed_bin(press_bin<=Zeu));
% 
% newtype=fittype('a*x+b');
% m2=fit(x,y,newtype);
% Kd_inZeu=(-1)*m2.a;
% 
% % Kd_inZpd
% x = press_bin(press_bin<=Zpd);
% y = log(Ed_QCed_bin(press_bin<=Zpd));
% 
% newtype=fittype('a*x+b');
% m2=fit(x,y,newtype);
% Kd_inZpd=(-1)*m2.a;
% end

% Get variables
PRES = press_Ed;
irr_380 = [PRES, Ed];
irr_380 = array2table(irr_380, 'VariableNames', {'PRES', 'IRR_380'});
irr_380.qc_380 = flag_Ed;

i_bon_380 = find(irr_380.qc_380 == 1 | irr_380.qc_380 == 2);
NewIRR_380 = irr_380(i_bon_380, :);

% Polynomial fit to remove the noise
p = polyfit(NewIRR_380.PRES, log(NewIRR_380.IRR_380), 4);
y_fit = polyval(p, NewIRR_380.PRES);
residuals = log(NewIRR_380.IRR_380) - y_fit;
mean_380 = mean(residuals, 'omitnan');
sd_380 = std(residuals, 'omitnan');
lim_sd2_380 = 2 * sd_380;

flag3_380 = residuals < (mean_380 - lim_sd2_380) | residuals > (mean_380 + lim_sd2_380);
no_cloudy_380 = find(~flag3_380);
new380 = NewIRR_380(no_cloudy_380, :);

% Define depth of safety
ilim_380 = new380.IRR_380(end);
up_380 = find(new380.IRR_380 >= ilim_380 * 1);
newdata_380 = new380(up_380, :);

% Surface irradiance calculation
% first calculation from the shallowest measurements
lim_z380_up_1 = find(newdata_380.IRR_380 >= newdata_380.IRR_380(1) / exp(1));% values within the penetration depth (1/e)
p_z1 = polyfit(newdata_380.PRES(lim_z380_up_1), log(newdata_380.IRR_380(lim_z380_up_1)), 2);
surf_1_380 = exp(p_z1(3));

% Second calculation using surf_1
New3802 = [[0, surf_1_380, NaN]; table2array(newdata_380)];
lim_z380_up_2 = find(New3802(:, 2) >= New3802(1, 2) / exp(1));
p_z2 = polyfit(New3802(lim_z380_up_2, 1), log(New3802(lim_z380_up_2, 2)), 2);
surf_2_380 = exp(p_z2(3));

% Third calculation using surf_2
New3803 = [[0, surf_2_380, NaN]; table2array(newdata_380)];
lim_z380_up_3 = find(New3803(:, 2) >= New3803(1, 2) / exp(1));
p_z3 = polyfit(New3803(lim_z380_up_3, 1), log(New3803(lim_z380_up_3, 2)), 2);
surf_3_380 = exp(p_z3(3));

New380 = [[0, surf_3_380, NaN]; table2array(newdata_380)];

% Bin by 1 meter
PRES = New380(:, 1);
VALUE = New380(:, 2);
brk_pts = -0.5 : 1 : round(max(PRES, [], 'omitnan')) + 0.5;
[~, ~, grp] = histcounts(PRES, brk_pts);
result = [];
a = unique(grp);
for ii = 1:length(a)
    idx = (grp == a(ii));
    % 如果该层没有数据，跳过
    if sum(idx) == 0
        continue;
    end
    % 当前层的中位数（对 VALUE 求）
    medVal = median(VALUE(idx), 'omitnan');
    % 记录当前 bin（以 grp 编号为深度）和中位值
    result = [result; a(ii), medVal];
end

% for ii = unique_groups(1):unique_groups(end)
%     for j = 2:size(New380, 2)
%         idx = New380_grp == ii;
%         New380(idx, j) = median(New380(idx, j), 'omitnan');
%     end
% end
% 
% New380 = unique(New380, 'rows');
% New380(:, 1) = New380_grp;
New380 = result;

% Kd within zpd calculation
i_z380 = find(New380(:, 1) <= Zpd);
% Zm1_380 = mean(New380(i_z380, 1));

if length(New380(i_z380, 2)) >= 3 && round(New380(1, 2), 2) > round(New380(2, 2), 2)
    mean_out = mean(New380(i_z380, 2), 'omitnan');
    sd_out = 3 * std(New380(i_z380, 2), 'omitnan');
    i_in = i_z380(New380(i_z380, 2) <= (mean_out + sd_out) & New380(i_z380, 2) >= (mean_out - sd_out));
    lm1 = polyfit(New380(i_in, 1), log(New380(i_in, 2)), 1);
    
    if (1 - var(log(New380(i_in, 2)) - polyval(lm1, New380(i_in, 1))) / var(log(New380(i_in, 2)))) >= 0.90
        Kd_zpd = -lm1(1);
        SE_zpd = std(log(New380(i_in, 2)) - polyval(lm1, New380(i_in, 1))) / sqrt(length(i_in));
    else
        Kd_zpd = NaN;
        SE_zpd = NaN;
    end
else
    Kd_zpd = NaN;
    SE_zpd = NaN;
end


% compare with kw
kw = pop_kw(lambda);
if Kd_zpd <= kw
    Kd_zpd = nan;
    SE_zpd = nan;
end

% Compare with the global range:
% from Organelli et al., (2017), ESSD
if Kd_zpd<=0.010 || Kd_zpd>=0.546
    Kd_zpd = nan;
    SE_zpd = nan;
end
end