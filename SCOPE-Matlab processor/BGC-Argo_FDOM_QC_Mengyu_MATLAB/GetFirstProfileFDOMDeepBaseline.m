function [FirstProfileFDOMDeepBaseline,DepthOfBaseline] = GetFirstProfileFDOMDeepBaseline(RawFDOM3arrays,PressOfFDOM3arrays)
%
%    GetFirstProfileFDOMDeepBaseline detail information:
%
%    SCOPE 07/02/2025
%
%    for SCOPE pj FDOM QC
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    Method via: Organelli et al., 2017 ESSD, DOI: 10.5194/essd-9-861-2017
%    This function is based on data from the first three cycles of the buoy to obtain a 
%    950-1000 deep baseline, in preparation for the subsequent FDOM QC, according to 
%    Organelli et al., 2017 ESSD, Section 2.2;
% 
%    Input:
%    RawFDOM3arrays: Raw FDOM 3 profiles into one array, nx1, ppb
%    PressOfFDOM3arrays: The press of raw FDOM profiles, should be the same size like
%               FDOM, nx1, m
%
%    Output:
%    FirstProfileFDOMDeepBaseline: Deep baseline of FDOM from 950-1000m if
%       there is, mx1, ppb, if not, then return an empty array, []
%    DepthOfBaseline: depth information of the baseline, mx1, m, 
%       if FirstProfileFDOMDeepBaseline is empty, then return an empty array
%
%


index = find(PressOfFDOM3arrays>=950 & PressOfFDOM3arrays<=1000);
if isempty(index) | length(index)<=3
    FirstProfileFDOMDeepBaseline = [];
    DepthOfBaseline = [];
    warning('There is no deep CDOM baseline')
    return
end

FDOM_first3_raw = RawFDOM3arrays;
FDOM_first3_raw_press = PressOfFDOM3arrays;

index_nan = find(isnan(FDOM_first3_raw));% take out nan
FDOM_first3_raw_press(index_nan) = [];
FDOM_first3_raw(index_nan) = [];

[FDOM_first3_press_sort,index] = sort(FDOM_first3_raw_press);% sort order according to press
FDOM_first3_sort = FDOM_first3_raw(index);

% remove the outliners
FDOM_first3_sort_rmoutliers = FDOM_first3_sort;
[~,TF] = rmoutliers(FDOM_first3_sort,'movmedian',5);% outliers of the raw profile were removed
FDOM_first3_sort_rmoutliers(TF) = nan;

FirstProfileFDOMDeepBaseline = FDOM_first3_sort_rmoutliers(FDOM_first3_press_sort>=950 & ...
    FDOM_first3_press_sort<=1000);% this is the FDOM baseline 
DepthOfBaseline = FDOM_first3_press_sort(FDOM_first3_press_sort>=950 & ...
    FDOM_first3_press_sort<=1000);% this is the baseline press

% plot(FDOM_first3_sort,FDOM_first3_press_sort);hold on;
% plot(FDOM_first3_sort_rmoutliers,FDOM_first3_press_sort)

end

