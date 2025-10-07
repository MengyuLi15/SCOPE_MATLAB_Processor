function Zp = CalculationZpFromEd(Ed,press_Ed)
%
%    CalculationZpFromEd detail information:
%
%    SCOPE 24/2/2025
%    in MATLAB R2024b
%
%    for SCOPE pj
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    This function is for calculate the penetration depth of Ed412
%
%    Inputs：
%    Ed：Ed profile, Nx1, par profile, and don t have to QC in advance;
%    press_Ed：Press of Ed profile, Nx1
%
%    Output:
%    Zp: penetration depth of Ed412, m
%       

irr_temp = Ed;
press_temp = press_Ed;
press_temp(isnan(irr_temp)) = [];
irr_temp(isnan(irr_temp)) = [];

 if isempty(irr_temp)
    Zp = nan;
    return;
 end

%% Shapiro-Wilk DARK IDENTIFICATION 
% Check the data and perform the Shapiro-Wilk normality test
% Initializes the Shapiro-Wilk normality test result vector

ShapiroWilkforsRadiometry = nan(length(irr_temp) - 4, 1);

for l = 1:length(ShapiroWilkforsRadiometry)
    try
        % Extracts the current subarray and performs the Lilliefors test
        currentData = irr_temp(l:end);
        [H, pValue,SWstatistic] = swtest(currentData,0.05); % Lilliefors 
        ShapiroWilkforsRadiometry(l) = pValue; % save p
    catch
        % If there is an error, copy the previous value
        if l > 1
            ShapiroWilkforsRadiometry(l) = ShapiroWilkforsRadiometry(l - 1);
        else
            ShapiroWilkforsRadiometry(l) = NaN; % If it is the first value, set it to NaN
        end
    end
end
% put the flag 3 according to the Lilliefors test results
index_flag_temp = find(abs(ShapiroWilkforsRadiometry)>=0.05); % have a little bit change than the original one
if ~isempty(index_flag_temp)
    flag_temp_firstdepth = index_flag_temp(1,1)+4;
    Zp = press_temp(flag_temp_firstdepth);
else
    Zp = nan;
    return;
end

end
