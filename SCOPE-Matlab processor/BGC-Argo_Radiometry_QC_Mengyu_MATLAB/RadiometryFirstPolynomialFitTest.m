function flag_Radiometry = RadiometryFirstPolynomialFitTest(Radiometry,press_Radiometry,flag_Radiometry)
%
%    RadiometryFirstPolynomialFitTest detail information:
%
%    SCOPE 17/1/2025
%
%    for SCOPE pj PAR AND IRR QC 
%
%    MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%    Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Organelli et al., 2016 JTECH, DOI: 10.1175/JTECH-D-15-0193.1
%    This function is for the first 4-Polynomial Fit Test for irr and par
%    profiles, according to Organelli et al., 2016 JTECH, Figure 6, Step 2;
%
%    QC=1 good data
%    QC=2 probably good
%    QC=3 probably bad 
%    QC=4 bad
%
%    Inputs：
%    Radiometry：Radiometry profile, Nx1, irr or par profile, and don t
%       have to remove flag 3 in advance;
%    press_Radiometry：Press of Radiometry profile, Nx1
%    flag_Radiometry：Initialize the flag before the first 4-Polynomial Fit
%       Test, Nx1
%
%    Output:
%    flag_Radiometry: flag after the first 4-Polynomial Fit Test, Nx1


Radiometry_forfit = Radiometry;
press_Radiometry_forfit = press_Radiometry;
index_all = (1:1:length(Radiometry))';

flag_3 = find(flag_Radiometry==3|flag_Radiometry==4);
flag_not3 = find(flag_Radiometry~=3 & flag_Radiometry~=4);

Radiometry_forfit(flag_3) = [];
press_Radiometry_forfit(flag_3) = [];
index_all(flag_3) = [];

if isempty(Radiometry_forfit)
    flag_Radiometry(isnan(flag_Radiometry)) = 4;
    return;
end

x = press_Radiometry_forfit;
y = log(Radiometry_forfit);
[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly4' ); % fittype: fourthdegree polynomial
[~, gof,output] = fit( xData, yData, ft );
R2_Poly4Overall = gof.rsquare;

if R2_Poly4Overall < 0.995
    flag_Radiometry(isnan(flag_Radiometry)) = 3;
else
    fit_residuals = output.residuals;
    mean_fit_residuals = mean(fit_residuals,'all',"omitnan");
    std_fit_residuals = std(fit_residuals,0,"all");
    flag_3_limite = 2.*std_fit_residuals;
    first_outliners_forfit = find(fit_residuals > (mean_fit_residuals+flag_3_limite) |...
        fit_residuals < (mean_fit_residuals-flag_3_limite));
    if ~isempty(first_outliners_forfit)
        first_outliners = flag_not3(first_outliners_forfit);
        flag_Radiometry(first_outliners) = 3;
    end
end

end