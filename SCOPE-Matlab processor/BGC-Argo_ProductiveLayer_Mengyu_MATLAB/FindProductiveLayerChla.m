function [Zpl, residual] = FindProductiveLayerChla(press,chla)
% FindProductiveLayerChla: A summary of this function is shown here
%       To find the bottom of the productive layer according to Rasse and
%       Dall'Olmo, 2019.
%       The defination of the bottom of the productive layer is 0.05 mg m-3 of Chla:
%           "We used chl to define the productive layer where living phytoplankton are 
%           present because our technique to compute export requires that no production 
%           is taking place below a given layer of the water column. We found that zp was 
%           always deeper than the MLD with average values (±standard deviation, WMO number)
%           of 139m (±12, 6901174) and 130m (±13, 6901175); zp approximately followed the 
%           26.6 kgm−3 isopycnal (Figures 2 and 3)."

    index = isnan(press) | isnan(chla);
    press(index) = []; chla(index) = [];

    if isempty(chla)
        Zpl = nan;
        residual = nan;
        return;
    end
    
    % remove noise
    [~,index] = rmoutliers(chla,"movmedian",11);
    press(index) = []; chla(index) = [];

    [~,b1] = sort(abs(chla - 0.05));% b1 = b1(1:50);
    [~,b2] = max(chla);
    index = b1 > b2;
    b = b1(index);

    if isempty(b)
        Zpl = nan;
        residual = nan;
        return;
    end

    Zpl = press(b(1));
    residual = chla(b(1)) - 0.05;
    
    % plot(press,chla)
end