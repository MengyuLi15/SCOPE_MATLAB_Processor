%    Readme
%
%    SCOPE 21/1/2025
%    in MATLAB R2024b
%
%    for SCOPE pj
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%
% This folder contains MATLAB QC programs for BGC-Argo's Ed and PAR.
% The main reference is Organelli et al., 2016, JTECH, DOI: 10.1175/ Jtech-D-15-0193.1
% 
% The main program is Radiometry_QC.m file
% The rest are necessary functions
%
%
% The processes include:
%
% 1. Read data
% 2. Preliminary QC
% 3. Fine correction
% 4. Calculate Zeu and Zpd based on PAR (or surface Chla estimates if no PAR is available)
% 5. Kd was calculated according to Ed and Zeu
