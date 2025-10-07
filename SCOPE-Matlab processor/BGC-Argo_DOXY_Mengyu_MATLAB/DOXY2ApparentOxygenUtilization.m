function AOU_profile = DOXY2ApparentOxygenUtilization(S_profile, T_profile, O2_meas_profile)
%
%    DOXY2ApparentOxygenUtilization detail information:
%
%    SCOPE 12/2/2025
%
%    for SCOPE pj: FDOM and DOXY
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Weiss, 1970, doi: https://doi.org/10.1016/0011-7471(70)90037-9
%    This function is for calculation of AOU (Apparent Oxygen Utilization) for 
%       the entire profile for BGC-Argo
%
% 
%    Input:
%    S_profile       - Salinity profile (PSS-78), nx1
%    T_profile       - Temperature profile (°C), nx1
%    O2_meas_profile - Dissolved oxygen profile (μmol/kg), nx1
%
%    Output：
%    AOU_profile     - AOU profile (μmol/kg), nx1



%    Make sure the input data is in the same shape
if ~isequal(size(S_profile), size(T_profile), size(O2_meas_profile))
    error('The dimensions of the input data must match!');
end

% Convert temperature (°C) to Kelvin (K)
T_K = T_profile + 273.15;

% Weiss (1970) formula coefficients
A1 = -173.4292; A2 = 249.6339; A3 = 143.3483; A4 = -21.8492;
B1 = -0.033096; B2 = 0.014259; B3 = -0.0017000;

% Calculated Oxygen saturation (μmol/kg)
ln_O2_sat = A1 + A2 .* (100 ./ T_K) + A3 .* log(T_K / 100) + A4 .* (T_K / 100) + ...
            S_profile .* (B1 + B2 .* (T_K / 100) + B3 .* (T_K / 100).^2);
O2_sat_profile = exp(ln_O2_sat) * 44.6596;  % Convert mL/L to μmol/kg

% Calculate the AOU profile
AOU_profile = O2_sat_profile - O2_meas_profile;
end