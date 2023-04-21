function semi_norm = penalisation_cost(Fx,norm_type)



% compute the TV amount
if strcmpi(norm_type, 'L1')

    % compute the element-wise norms
    Fx = abs(Fx);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'L2')

    % compute the block-wise norms
    Fx = Fx.^2;
    Fx = sum(Fx, 1);
    Fx = sqrt(Fx);
    Fx = shiftdim(Fx,1);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'inf')
    
    % compute the block-wise norms
    Fx = abs(Fx);
    Fx = max(Fx, [], 1);
    Fx = shiftdim(Fx,1);
    
    % sum the norms
    semi_norm = sum( Fx(:) );

elseif strcmpi(norm_type, 'nuclear')
    
    % compute the block-wise norms
    [u s v] = sv_dec(Fx);
    Fx = abs(s);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'frobenius')
   
    % compute the block-wise norms
    Fx = Fx.^2;
    Fx = sum(Fx, 1);
    Fx = sum(Fx, 2);
    Fx = sqrt(Fx);
    Fx = squeeze(Fx);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
elseif strcmpi(norm_type, 'spectral')
    
    % compute the block-wise norm
    [u s v] = sv_dec(Fx);
    Fx = max(s, [], 1);
    Fx = squeeze(Fx);
    
    % sum the norms
    semi_norm = sum( Fx(:) );
    
else
    error('Norm not recognized.') 
end
