function [Variable_QCed,Variable_flag] = BGCArgo_ncread_VariableAndADJUSTED(filename,bgc_variables,VariableName,VariableNameQC,VariableName_ADJUSTED,VariableName_ADJUSTED_QC)
%
%    BGCArgo_ncread_VariableAndADJUSTED detail information:
%
%    SCOPE 17/1/2025
%
%    for SCOPE pj FDOM
%
%    MENGYU LI, CNR-ISMAR, mengyu.li@artov.ismar.cnr.it
%    Emanuele Organelli, PI, emanuele.organelli@cnr.it
%
%    For a BGC float, in the example of Chla, the .nc file is likely to have two variable forms
%    One is the original form, like Chla, and the other may be variables that have been adjusted 
%      in the data center, like Chla_ADJUSTED
%    The original form has not been processed in the data center and needs to focus on QC analysis
%    This function by input nc file info data and required variables, read the corresponding 
%      data and QC file, do the preliminary QC
%
% 
%    Input:
%    bgc_variables: info of the nc file read earlier, generally recommended
%                   bgc_variables = {ncinfo(filename).Variables.Name}';
%    VariableName: The target variable name, as a string, like 'PRES'
%    VariableNameQC: The target variable name corresponds to the QC flag, as a string, 
%                    like 'PRES_QC'
%    VariableName_ADJUSTED: A variable that has been adjusted for the desired target, 
%                           as a string, like 'PRES_ADJUSTED'
%    VariableName_ADJUSTED_QC: The corresponding QC flag of the adjusted variable 
%                              that is required, as a string, similar to 'PRES_ADJUSTED_QC'
%
%    Output:
%    Variable_QCed: QC data of the target variable, consistent with the size read in; 
%                   If this variable is not present and Variable_flag=0 is returned, 
%                   the variable returns empty
%    Variable_flag: Returns whether a required variable is available, 
%                   with 1 indicating that the variable is present and 
%                   0 indicating that the variable is not present

flag_orginal = any(strcmp(bgc_variables,VariableName));
flag_ADJUSTED = any(strcmp(bgc_variables,VariableName_ADJUSTED)); 
    if flag_ADJUSTED % we have ADJUSTED press, that could be the best,then use the adjusted
        Variable_QCed = ncread(filename,VariableName_ADJUSTED);
        Variable_QC = ncread(filename,VariableName_ADJUSTED_QC);
        Variable_QCed(Variable_QC==3|Variable_QC==4|isnan(Variable_QC)|Variable_QC==6|...
            Variable_QC==7|Variable_QC==9) = NaN;
        Variable_flag = true;
    elseif flag_orginal
        Variable_QCed = ncread(filename,VariableName);
        Variable_QC = ncread(filename,VariableNameQC);
        Variable_QCed(Variable_QC==3|Variable_QC==4|isnan(Variable_QC)|Variable_QC==6|...
                Variable_QC==7|Variable_QC==9) = NaN;
    
        Variable_flag = true;
    else
        Variable_QCed = [];
        Variable_flag = false;
    end
end