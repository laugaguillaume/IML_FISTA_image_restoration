function [R,P] = get_information_transfer(parameters,level)
%GET_INFORMATION_TRANSFER Summary of this function goes here
R_index = matlab.lang.makeValidName(['level' num2str(level-1) 'R']); P_index = matlab.lang.makeValidName(['level' num2str(level-1) 'P']);
R = parameters.blur.(R_index); P = parameters.blur.(P_index);
end

