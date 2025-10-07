function [flag_Radiometry, isLessThan5Points] = RadiometryProfilePointsNumTest(Radiometry,press_Radiometry,flag_Radiometry)
%
%    RadiometryProfilePointsNumTest detail information:
%
%    SCOPE 17/1/2025
%    in MATLAB R2024b
%
%    for SCOPE pj PAR AND IRR QC 
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Organelli et al., 2016 JTECH, DOI: 10.1175/JTECH-D-15-0193.1
%    This function is for the Profile point number test before the first 4-Polynomial Fit 
%    Test for irr and par profiles, according to Organelli et al., 2016 JTECH, 
%    Figure 6, Step 1;
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
%    flag_Radiometry: flag after the fprofile points number Test, Nx1
%    isLessThan5Points: flag=1 means the profile points number less 5,
%       flag=0 means the profile number is more than 5.

    Radiometry_fortest = Radiometry;
    press_Radiometry_fortest = press_Radiometry;
    flag_34 = find(flag_Radiometry==3|flag_Radiometry==4);
    Radiometry_fortest(flag_34) = [];
    press_Radiometry_fortest(flag_34) = [];

    if length(Radiometry_fortest) > 5
        isLessThan5Points = 0;
    else
        flag_Radiometry(isnan(flag_Radiometry))=3;
        isLessThan5Points = 1;
    end

end