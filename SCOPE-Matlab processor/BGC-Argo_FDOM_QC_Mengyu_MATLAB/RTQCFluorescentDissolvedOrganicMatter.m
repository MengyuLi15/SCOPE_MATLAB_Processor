function [FDOM_adjusted,Press_Raw,FDOM_QCFlag] = RTQCFluorescentDissolvedOrganicMatter(FDOM_Raw,Press_Raw,FirstProfileFDOMDeepBaseline)
%
%    RTQCFluorescentDissolvedOrganicMatter detail information:
%
%    SCOPE 27/1/2025
%
%    for SCOPE pj FDOM QC
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Organelli et al., 2017 ESSD, DOI: 10.5194/essd-9-861-2017
%    This function performs QC on FDOM from BGC-Argo
%    profiles, according to Organelli et al., 2017 ESSD, Section 2.2;
%
% 
%    Input:
%    FDOM_Raw: FDOM raw profiles without QC, nx1, ppb
%    Press_Raw: The press of FDOM profile, should be the same size like
%               FDOM, nx1, m
%    FirstProfileFDOMDeepBaseline: We believe that CDOM and salinity are conserved over 
%               950-1000m, if it is the same water mass. 
%               Therefore, the deep CDOM of the first three profiles was collected as the 
%               baseline for drift of the profiles
%
%    Output:
%    FDOM_adjusted: FDOM data after QC, including the result of float drift correction, nx1, ppb
%    Press_Raw: Exactly the same depth data, in order to remind the user of depth，nx1, m
%    FDOM_QCFlag：Quality control flag
%
%    QC=1 good data
%    QC=2 probably good
%    QC=3 probably bad 
%    QC=4 bad
% for test
% FDOM_Raw = FDOM_temp;
% Press_Raw = Press_temp;
% FirstProfileFDOMDeepBaseline = FirstProfileFDOMDeepBaseline;


    Press_temp = Press_Raw; % initialize
    FDOM_temp = FDOM_Raw;
    FDOM_QCFlag = nan(length(FDOM_temp),1);
    index_all = (1:1:length(FDOM_temp))';
%----------------0 step:  if there is no deep baseline, then we cannot corrent drift------
    FirstProfileFDOMDeepBaseline(isnan(FirstProfileFDOMDeepBaseline)) = [];
    if isempty(FirstProfileFDOMDeepBaseline)
        FDOM_QCFlag(isnan(FDOM_QCFlag)) = 3;
        FDOM_adjusted = FDOM_Raw;
        warning('There is no deep CDOM baseline in this wmo, flag = 3')
        return
    end

    index_deep = Press_temp>=950 & Press_temp<=1000;  
    deepFDOM = FDOM_temp(index_deep);
    deepFDOM(isnan(deepFDOM)) = [];
    if isempty(deepFDOM)
        FDOM_QCFlag(isnan(FDOM_QCFlag)) = 3;
        FDOM_adjusted = FDOM_Raw;
        warning('There is no deep CDOM baseline in this profile, flag = 3')
        return
    end

%----------------First step:  specific range test from Wetlab user manual and nan test----
%-----via: http://gyre.umeoce.maine.edu/resources/instrumentation/ECO_sensors_index.pdf---
%     FDOM 370/460 nm Range: 0-375, Sensitivity: 0.184 ppb https://doi.org/10.5194/essd-9-861-2017
    index_outnan = (FDOM_temp>5 | FDOM_temp < 0) | isnan(FDOM_temp) ;
    FDOM_QCFlag(index_outnan) = 4;

    Press_temp(index_outnan) = [];
    FDOM_temp(index_outnan) = [];
    index_all(index_outnan) = [];
%-----------------Second step: spikes outside the 25th and 75th quantiles test------------
    y = smoothdata(FDOM_temp,'movmedian','SmoothingFactor',0.25);
    spikes = FDOM_temp-y;
    % histogram(spikes)
    index_2575 = find(spikes<quantile(spikes,0.25)|spikes>quantile(spikes,0.75));
    % [~,index_2575] = rmoutliers(FDOM_temp,"percentiles",[25 75]);
    FDOM_QCFlag(index_all(index_2575)) = 4;

    Press_temp(index_2575) = [];
    FDOM_temp(index_2575) = [];
    index_all(index_2575) = [];

%-----------------Third step: absolute residual >4 calculated with mean filter -----------
    [~,index_residual] = rmoutliers(FDOM_temp,"movmedian",[3 3]);% hampel filter
    FDOM_QCFlag(index_all(index_residual)) = 4;

    Press_temp(index_residual) = [];
    FDOM_temp(index_residual) = [];
    index_all(index_residual) = [];
    % plot(FDOM_temp,Press_temp)
%-------Forth step: a smooth filter first a median filer 5 and then an average filter 7---
    FDOM_Raw(FDOM_QCFlag==4|FDOM_QCFlag==3) = nan;
    FDOM_Raw = smoothdata(FDOM_Raw,"movmedian",5,"omitnan");
    FDOM_Raw = smoothdata(FDOM_Raw,"movmean",7,"omitnan");

%------------------------Fifth step: deep baseline drift correction-----------------------
    index_deep = find(Press_temp>=950 & Press_temp<=1000);
    if isempty(index_deep)
        FDOM_QCFlag(isnan(FDOM_QCFlag)) = 3;
        FDOM_adjusted = FDOM_Raw;
        warning('There is no deep CDOM in this profile, so flag = 3')
        return
    else
        deep_FDOM =  FDOM_temp(index_deep); 
        deep_FDOM_mean = mean(deep_FDOM,'all','omitnan');
        offset = deep_FDOM_mean-mean(FirstProfileFDOMDeepBaseline,'all','omitnan');
        FDOM_adjusted = FDOM_Raw-offset;
    end

    FDOM_QCFlag(isnan(FDOM_QCFlag)) = 1;

end

