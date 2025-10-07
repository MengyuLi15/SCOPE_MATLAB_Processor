function [var_QCED_1d,press_var_QCED_1d] = BGCArgo_single_prof_QC(press_QCED_2d,var_raw_2d,var_QC_2d)
%
%    BGCArgo_single_prof_QC detail information:
%
%    Programed by MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%                 Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%
%    This function is for the RT Quality Control for irr and par BGC-Argo
%
%   BGCArgo_single_prof_QC，Used to conduct preliminary flag QC for the profile of a 
%   single variable for a single float file and generate 1d depth-dependent profiles and depths
%
%   [var_QCED_1d,press_var_QCED_1d] = BGCArgo_single_prof_QC(press_QCED_2d,var_raw_2d,var_QC_2d)
%
%   Inputs：
%   press_QCED_2d：Press after QC in 2D
%   var_raw_2d：The raw data waiting of QC, 2D, Read from nc
%   var_QC_2d：The raw data QC file, 2D, Read from nc
%
%    Outputs：
%    var_QCED_1d：Data passing through QC with increasing depth, 1D
%    press_var_QCED_1d：Corresponding depth of the data, 1D

%     PSAL = ncread(filename_R,'PSAL'); read the var_raw_2d
    size_mat = size(var_raw_2d); % get the var 2d row and lines
    var_QCED_2d = NaN(size_mat(1),size_mat(2)); % generate the empty 2d qced var
%     PSAL_QC = ncread(filename_R,'PSAL_QC');  read the QC flags var_QC_2d

    for roll = 1:size_mat(2)
        var_accept = var_raw_2d(:,roll);
        flag = double(string(var_QC_2d(:,roll)));
        index_accept = find(flag==3|flag==4|flag==9);
        if isempty(index_accept)
            var_QCED_2d(:,roll) = var_accept;
            continue;
        else
            var_accept(index_accept) = NaN; 
            var_QCED_2d(:,roll) = var_accept;
        end
    end
    % 这段是为了找到对应深度的那几列 这会得到的也潜在是2d的剖面 或者最好是已经是1d的了
    index_sal = find(~all(isnan(var_QC_2d))==true);
    PRES_QCED_var = press_QCED_2d(:,index_sal);
    s = size(PRES_QCED_var);
    var_QCED_2d = var_QCED_2d(:,index_sal);
    % 一般有意义都得是正的 不是的话给0 给nan也行 看情况这里先给的0
    var_QCED_2d(var_QCED_2d<0) = nan;
    % 有可能有几列对应的
    PRES_QCED_var = reshape(PRES_QCED_var,[s(1)*s(2) 1]);
    var_QCED_2d = reshape(var_QCED_2d,[s(1)*s(2) 1]);
    % 现在都变成一列的了 包括var和深度
    index = isnan(var_QCED_2d)|isnan(PRES_QCED_var); % 去掉所有的nan值
    var_QCED_2d(index) = [];PRES_QCED_var(index) = [];
    % 再去一遍nan
    index = isnan(PRES_QCED_var);var_QCED_2d(index) = [];PRES_QCED_var(index) = [];
    % 传感器可能上下都测了 现在要根据深度排序
    [PRES_QCED_var,index] = sort(PRES_QCED_var,'ascend');
    var_QCED_2d = var_QCED_2d(index);
    [PRES_QCED_var,index] = unique(PRES_QCED_var);
    var_QCED_2d = var_QCED_2d(index,:);
    var_QCED_1d = var_QCED_2d;
    press_var_QCED_1d = PRES_QCED_var;


end

