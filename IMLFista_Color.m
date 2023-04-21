function [X,X_iteration,fun,counter_descent,coarse_descent_f,time_steps,snr_steps]=IMLFista_Color(level,fh,gh,gradfh,proxgh,X,Xref,parameters_reg,parameters,lmax,bk_fine,bk_coarse,acceleration) 
%% V-cycle ?
% Inputs : 
% Outputs : 
% x : image
% iteration : iterations number
% fun : objective function values at each iteration
% grad : gradient norm values at each iteration
% counter_descent : number of iterations where coarse models were called
% time_steps : CPU time steps
% X_iterations : evolution of X at given iterations
if nargin < 13
    acceleration = 'ISTA';
end
if level>1
    it_number = matlab.lang.makeValidName(['maxit_level' num2str(level)]); 
    max_iteration  = parameters.(it_number);
    % Definition of variables we gather to plot
    fun           = NaN(2,max_iteration); 
    time_steps    = NaN(1,max_iteration);
    time_steps(1) = 0;
    snr_steps = NaN(1,max_iteration);
    snr_steps(1) = snrdeblur(X,Xref);
    coarse_descent_f = NaN(1,max_iteration);
    % Extract functions, gradients, proximal operator handles and
    % information transfer matrix
    [fH,gH,gradfH,proxgH,DH,DH_adjoint] = get_func_FB_D(parameters,level);
    proxg = parameters.proximal.proxg;
    [R,P] = get_information_transfer(parameters,level);
    %
    alpha_index = matlab.lang.makeValidName(['level' num2str(level) 'alpha']);
    alpha_lip   = parameters.(alpha_index);
    %
    Fh          = @(U) fh(U) + gh(U);
    [N,M]       = size(X);
    NH = N/2;MH = M/2;
    % Initialization of the descent
    iteration      = 0;
    fun(:,1)         = [fh(X);gh(X)];
    gradfx         = gradfh(X);
    stop_criterion = iteration<max_iteration;
    Y              = X;
    %last_coarse_iteration = zeros(N,M);
    counter_descent    = 0;
    got_down_once      = 0;
    enter_coarse_model = 1;
    saved_coarse = 0;
    while stop_criterion
        tic
        Xm = X;
        if got_down_once>= parameters.p %|| (level<lmax && iteration < 1 )  
            enter_coarse_model = 0;
        else
            X = Y;
            if iteration == 0
                DhXh = parameters.dir_op(X);
                proxgk = proxg(DhXh,alpha_lip);
                direction_h = (1/alpha_lip).*parameters.adj_op(DhXh-proxgk)+gradfx;
            else
                DhXh = parameters.dir_op(X);
                proxgk = proxg(DhXh,alpha);
                direction_h = (1/alpha).*parameters.adj_op(DhXh-proxgk)+gradfx;
            end
            directionRh = D3RP(direction_h,R);
        end
        if enter_coarse_model % we go down one level
            counter_descent = counter_descent +1; % increment each time we go down 
            Xm = X;
            xH0 = D3RP(X,R);  
            DHxH0 = DH(xH0);
            vH = directionRh - gradfH(xH0)-parameters.moreau^(-1).*DH_adjoint(DHxH0-proxgH(DHxH0,parameters.moreau));
            mH = @(x_coarse) fH(x_coarse) + vH(:)'*(x_coarse(:)-xH0(:));
            gradmH = @(x_coarse) gradfH(x_coarse) + vH;
            % Coarse model
            [xH] = Coarse_Model_GMoreauColor(level-1,mH,gH,gradmH,proxgH,DH,DH_adjoint,xH0,parameters,lmax,bk_coarse);
            % Coarse correction term
            step_size = 1;
            s = D3RP(xH-xH0,P); 
            while Fh(Y + step_size.*s)>Fh(Y)
                step_size = step_size/2;
            end
            Y = Y + step_size.*s;
            gradfx = gradfh(Y);
            if level<=lmax
                got_down_once = got_down_once+1;
            end
        end
        [alpha,ind,X] = backtracking_ForwardBackward_Color(alpha_lip,Y,-gradfx,gradfx,fh,gradfh,proxgh,bk_fine,parameters_reg);
        beta = (iteration+1)/(iteration+1+4);
        Y = X + beta*(X-Xm);
        iteration = iteration + 1;
        gradfx = gradfh(Y);
        fun(:,iteration+1)=[fh(X);gh(X)];
        if sum(fun(:,iteration))<sum(fun(:,iteration+1))
            parameters_reg.epsilon = parameters_reg.epsilon/10;
        end
        if enter_coarse_model && ind
            coarse_descent_f(iteration+1) = fun(iteration+1);
        end
        t = toc;
        stop_criterion = iteration<max_iteration;
        time_steps(iteration+1) = t + time_steps(iteration);
        snr_steps(iteration+1) = snrdeblur(X,Xref);
        if got_down_once >= parameters.p && saved_coarse<1
            X_iteration = X;
            saved_coarse = saved_coarse +1;
        end
    end

else % called if we have only one level
    switch acceleration
        case 'ISTA'
            [X,Xk,fun,time_steps,snr_steps]=StandardForwardbackward_Color(fh,gh,gradfh,proxgh,parameters_reg,parameters,X,Xref,level,bk_fine);
            X_iteration = Xk; counter_descent = 0; coarse_descent_f = 0;
        case 'FISTA'
            [X,Xk,fun,time_steps,snr_steps]=StandardFISTA_Color(fh,gh,gradfh,proxgh,parameters_reg,parameters,X,Xref,level,bk_fine);
            X_iteration = Xk; counter_descent = 0; coarse_descent_f = 0;
    end
end

end

