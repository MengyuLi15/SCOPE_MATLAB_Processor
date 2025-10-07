function [Zeu,Zpd] = CalculationZeuZpdFromPARProfile(PAR,press_PAR,flag_PAR)
%
%    CalculationZeuZpdFromPARProfile detail information:
%
%    SCOPE 20/1/2025
%    in MATLAB R2024b
%
%    for SCOPE pj
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    This function is for calculate the euphotic depth Zeu and first
%    optical depth Zpd from the PAR profile
%
%    Inputs：
%    PAR：PAR profile, Nx1, par profile, and don t
%       have to QC in advance;
%    press_PAR：Press of PAR profile, Nx1
%    flag_PAR：PAR QC flag, Nx1
%
%    Output:
%    Zeu: euphotic depth Zeu, m
%    Zpd: first optical depth Zpd, m
%       
% PAR_surface = mean(PAR(press_PAR<=3),"all",'omitnan');% surface PAR
% PAR_surface_1percent = PAR_surface.*0.01;% 1% PAR
% PAR = PARt;
% press_PAR = press_PARt;
% flag_PAR = ones(length(PARt),1);

PAR_QCed = PAR;
press_PAR_QCed = press_PAR;
PAR_QCed(flag_PAR==3|flag_PAR==4)  = nan;
press_PAR_QCed(isnan(PAR_QCed),:) = [];
PAR_QCed(isnan(PAR_QCed),:) = [];% after the qc and then for the next calculation

flag_PAR(flag_PAR==3 | flag_PAR==4) = [];

 if isempty(PAR_QCed)
    Zeu = nan;
    Zpd = nan;
    return;
 end

 %% 1.Polynomial fit to remove othe noise
irr_PAR = [press_PAR_QCed, PAR_QCed];
irr_PAR = array2table(irr_PAR, 'VariableNames', {'PRES', 'IRR_PAR'});
irr_PAR.qc_PAR = flag_PAR;
% Select valid QC values (1 or 2)
i_bon_PAR = find(irr_PAR.qc_PAR == 1 | irr_PAR.qc_PAR == 2);
NewPAR = irr_PAR(i_bon_PAR, :);

% Polynomial fit to remove noise
P = polyfit(NewPAR.PRES, log(NewPAR.IRR_PAR), 4);
fitted_values = polyval(P, NewPAR.PRES);
residuals = log(NewPAR.IRR_PAR) - fitted_values;

mean_PAR = mean(residuals, 'omitnan');
sd_PAR = std(residuals, 'omitnan');
lim_sd2_PAR = 2 * sd_PAR;

flag3_PAR = residuals < (mean_PAR - lim_sd2_PAR) | residuals > (mean_PAR + lim_sd2_PAR);% out of fit 
no_cloudy_PAR = ~flag3_PAR;
newPAR = NewPAR(no_cloudy_PAR, :);

% Set depth of safety
ilim_PAR = newPAR.IRR_PAR(end);
up_PAR = find(newPAR.IRR_PAR >= ilim_PAR * 1);
newdata_PAR = newPAR(up_PAR, :);

% First Optical Depth Calculation
lim_zpd_up_1 = find(newdata_PAR.IRR_PAR >= newdata_PAR.IRR_PAR(1) * 0.01);
lim_1 = newdata_PAR.PRES(end);
Zpd1 = lim_1 / 4.6;
i_lim1 = find(newdata_PAR.PRES <= Zpd1);

P1 = polyfit(newdata_PAR.PRES(i_lim1), log(newdata_PAR.IRR_PAR(i_lim1)), 4);
surf_1 = exp(P1(end));

% Second Calculation
NewPAR2 = [[0, surf_1, NaN]; table2array(newdata_PAR)];
lim_zpd_up_2 = find(NewPAR2(:,2) >= NewPAR2(1,2) * 0.01);
lim_2 = NewPAR2(lim_zpd_up_2(end), 1);
Zpd2 = lim_2 / 4.6;
i_lim2 = find(NewPAR2(:,1) <= Zpd2);
P2 = polyfit(NewPAR2(i_lim2,1), log(NewPAR2(i_lim2,2)), 2);
surf_2 = exp(P2(end));

% Third Calculation
NewPAR3 = [[0, surf_2, NaN]; table2array(newdata_PAR)];
lim_zpd_up_3 = find(NewPAR3(:,2) >= NewPAR3(1,2) * 0.01);
lim_3 = NewPAR3(lim_zpd_up_2(end), 1);
Zpd3 = lim_3 / 4.6;
i_lim3 = find(NewPAR3(:,1) <= Zpd3);
P3 = polyfit(NewPAR3(i_lim3,1), log(NewPAR3(i_lim3,2)), 2);
surf_3 = exp(P3(end));

% Final Optical Depth Calculation
NEWPAR = [[0, surf_3, NaN]; table2array(newdata_PAR)];
Residual = abs(NEWPAR(:,2) - NEWPAR(1,2) * 0.01);
[~,index] = min(Residual,[],'omitnan');
% lim_zpd_up_4 = find(NEWPAR(:,2) >= NEWPAR(1,2) * 0.01);
lim_4 = NEWPAR(index, 1);
Zpd = lim_4 / 4.6;

Zeu = lim_4;

if Zeu<10.5 || Zeu>203.8
    Zeu = nan; Zpd = nan;
end

end