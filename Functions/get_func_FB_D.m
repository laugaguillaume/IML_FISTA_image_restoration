function [fH,gH,gradfH,proxgH,DH,DH_adjoint] = get_func_FB_D(parameters,level)
%GET_FUNC_FB EXTRACT FUNCTIONS HANDLES FOR CODE READABILITY
%   Detailed explanation goes here
fH_index = matlab.lang.makeValidName(['level' num2str(level-1) 'fH']) ; gH_index = matlab.lang.makeValidName(['level' num2str(level-1) 'gH']);
gradfH_index = matlab.lang.makeValidName(['level' num2str(level-1) 'gradfH']);
%
fH = parameters.functions.(fH_index); gH = parameters.functions.(gH_index); 
gradfH = parameters.gradient.(gradfH_index);
%
proxgH_index = matlab.lang.makeValidName(['level' num2str(level-1) 'proxgH']);
proxgH = parameters.proximal.(proxgH_index);
%
DH_index = matlab.lang.makeValidName(['level' num2str(level-1) 'dir_opH']);
DH = parameters.op.(DH_index);
%
DH_adjoint_index = matlab.lang.makeValidName(['level' num2str(level-1) 'adj_opH']);
DH_adjoint = parameters.op.(DH_adjoint_index);
end

