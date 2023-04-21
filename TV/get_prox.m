function [prox] = get_prox(norm_type)
%GET_PROX Summary of this function goes here
%   Detailed explanation goes here
if strcmpi(norm_type, 'L1ball')
    
    prox = @(y,gamma) project_L1(y, 1, 1);
 
elseif strcmpi(norm_type, 'L1reg')
    
    prox = @(y,gamma) prox_abs(y, 1);
    
elseif strcmpi(norm_type, 'L2ball')
    
    prox = @(y,gamma) project_L1L2(y, 1, 1);

elseif strcmpi(norm_type, 'L2reg')
    
    prox = @(y,gamma) prox_L2(y, 1, 1);    
    
elseif strcmpi(norm_type, 'infball')
    
    prox = @(y,gamma) project_L1Linf(y, 1, 1);

elseif strcmpi(norm_type, 'nuclearball')
    
    prox = @(y,gamma) project_L1_Nuclear(y, 1, 1);
    
elseif strcmpi(norm_type, 'nuclearreg')
    
    prox = @(y,gamma) prox_L1_Nuclear(y, 1, 1);    
    
elseif strcmpi(norm_type, 'frobeniusball')
        
    prox = @(y,gamma) project_L1_Frobenius(y, 1, 1);
    
elseif strcmpi(norm_type, 'spectralball')
    
    prox = @(y,gamma) project_L1_Spectral(y, 1, 1);
end

end

