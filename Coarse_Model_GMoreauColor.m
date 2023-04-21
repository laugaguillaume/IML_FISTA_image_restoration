function [X]=Coarse_Model_GMoreauColor(level,f,g,gradf,proxg,D,D_adjoint,X,parameters,lmax,bk) 
%%
% Inputs : 
% Outputs : 
% x : image
% p : iterations number
% fun : objective function values at each iteration
% grad : gradient norm values at each iteration
if level>1
    it_number = matlab.lang.makeValidName(['maxit_level' num2str(level)]); 
    max_iteration = parameters.(it_number);
    % Definition of coarse model functions
    [fH,gH,gradfH,proxgH,DH,DH_adjoint] = get_func_FB_D(parameters,level);    
    [R,P] = get_information_transfer(parameters,level);
    %
    [N,M] = size(X);
    NH = N/2;MH = M/2;
    % Initialization of the descent
    iteration = 0;
    DX = D(X);
    proxgk = proxg(DX,parameters.moreau);
    Fx = f(X) + g(proxgk) + 1/(2*parameters.moreau)*norm(DX(:)-proxgk(:))^2;
    gradFx=gradf(X)+parameters.moreau^(-1).*D_adjoint(DX-proxgk);  
    stop_criterion = iteration<max_iteration;
    %last_coarse_iteration = zeros(N,M);
    counter_descent = 0;
    got_down_once = 0;
    alpha_lip = parameters.alpha_moreau;
    
    while stop_criterion
        % coarse level initialization (used if enter_coarse_model is
        % fulfilled
        if got_down_once>=1%parameters.p
            enter_coarse_model = 0;
        else
            gradFHx = D3RP(gradFx,R);
            enter_coarse_model=1;
            %enter_coarse_model= norm(reshape(X-last_coarse_iteration,N*M,1))>=parameters.eta*norm(reshape(last_coarse_iteration,N*M,1)) && norm(gradFHx)>parameters.eps ... 
            %    && norm(gradFHx)>=parameters.kappa*norm(gradFx(:)); % Condition 1
        end
        if enter_coarse_model % we go down one level
            xH0 = D3RP(X,R);
            DHxH0 = DH(xH0);
            counter_descent = counter_descent +1; % increment each time we go down 
            vH = gradFHx -gradfH(xH0)-parameters.moreau^(-1).*DH_adjoint(DHxH0-proxgH(DHxH0,parameters.moreau));
            mH = @(x_coarse) fH(x_coarse) + vH(:)'*(x_coarse(:)-xH0(:));
            gradmH = @(x_coarse) gradfH(x_coarse) + vH;
            % Coarse model
            xH = Coarse_Model_GMoreauColor(level-1,mH,gH,gradmH,proxgH,DH,DH_adjoint,xH0,parameters,lmax,bk);
            % Coarse correction term
            s = D3RP(xH-xH0,P);
            got_down_once = got_down_once + 1;
        else % usual gradient step
            s =-gradFx;
        end
        % Shape change to accomodate dot product
        % Search for optimal step size
        [alpha,ind,X,proxgk]=backtracking_GradientMoreauColor(alpha_lip,parameters.moreau,X,s,gradFx,Fx,f,g,proxg,D,D_adjoint,bk);
        DX = D(X);
        if bk ==1
            Fx = f(X) + g(proxgk) + 1/(2*parameters.moreau)*norm(DX-proxgk)^2;
        end
        gradFx=gradf(X)+parameters.moreau^(-1).*D_adjoint(DX-proxgk);
        iteration = iteration + 1;
        stop_criterion = iteration<max_iteration;
    end
elseif level == 1 && lmax>1
    X=Coarsest_GMoreauColor(f,g,gradf,proxg,D,D_adjoint,parameters,X,level,bk); 
    return
end
