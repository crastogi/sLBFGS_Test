function [path, function_values, grad_norm] = sLBFGS_k2SVRG(ff, x0, ...
    maxEpoch, stepsize, k, hessianBatch, hessianPeriod, memorySize, verbosity, x_star)
    global data nDataPoints;
    nDataPoints = 0;
    
    d = size(data, 2);
    N = size(data, 1);
    bH = hessianBatch;
    bs = 50;
    M = memorySize;
    eta = stepsize;
    epsilon = 1E-5;
    L = hessianPeriod;
    delta = 0;
    fdHVPStepSize = 5E-2;
    gradientNormBound = 100;
    % Additional k2-SVRG parameters
    l = ceil(N/(k*bs));
    S_l = l;
    
    % Number of currently computed Hessian correction pairs
    r = 0;
    currDepth = 0;
    hu_iters = 0;
    
    % Average of path travelled in the current and previous inverse hessian updates
    u_r = zeros(1, d);
    u_r_prev = zeros(1, d);
    % Components of two-loop update
    rho = zeros(1, M);
    s   = zeros(M, d);
    y   = zeros(M, d);
    % Store path, etc.
    path = [];
    function_values = [];
    grad_norm = [];
    
    % Need an size N array to store idx of locations. Initialize to 1
    theta_m_binding = ones(1, N);
    % Initialize theta matrix to 1 position set to x0
    theta_m = x0';
    % Set aBar_m
    [~, aBar_m] = ff(theta_m', 1:N);
    x_t = x0';
    xBar = x0';
    
    % Outer loop
    for m = 0:(maxEpoch-1)                             % Line 5
        % Print some information about the current epoch
        if (verbosity)
            disp(['Epoch: ' num2str(m) '; |aBar|: ' num2str(log10(norm(aBar_m))) '; Dist to *: ' num2str(log10(norm(xBar-x_star)))]);
        end
        
        ind = randperm(N);                      % Line 6
        bs_ind = randperm(N);
        for j = 0:(k-1)                         % Line 7
            x_t_pos = zeros(l,d);
            for t = 0:(l-1)                     % Line 9
                x_t_pos(t+1,:) = x_t;
                %i_t = randi(N);                 % Line 10
                %i_t = randsample(1:N, bs);      % Line 10
                i_t = bs_ind((j*l*bs+t*bs+1):(j*l*bs+(t+1)*bs));
                % Line 11
%                 [~, a_it] = ff(theta_m(theta_m_binding(i_t),:)', i_t);
%                 [~, f_it] = ff(x_t', i_t);
%                 g_it = f_it' - a_it' + aBar_m';
                sub_bindings = theta_m_binding(i_t);
                a_it = 0;
                for curr_theta = unique(sub_bindings)
                    [~, sg] = ff(theta_m(curr_theta,:)', i_t(sub_bindings==curr_theta));
                    a_it = a_it + sum(sub_bindings==curr_theta)*sg;
                end
                a_it = a_it/bs;
                [~, f_it] = ff(x_t', i_t);
                g_it = f_it' - a_it' + aBar_m';
                
                % Check to see if gradient vector has usable size
                if (norm(g_it)<1E-10)
                    % Potential convergence/scale issue; break
                    % TODO: fix
                    break
                end
                
                % Update u_r with current position
                u_r = u_r + x_t;
                % Compute next iteration step; condition the gradient so as
                % not to produce rapid oscilations @extreme function values
                if (r < 1)
                    % H_0 = I Until one hessian correction takes place
                    effGrad = g_it;
                else
                    % Compute the two-loop recursion product
                    effGrad = twoLoopRecursion(g_it);
                end
                %Bound the effective gradient update step
                egNorm = norm(effGrad);
                divisor = 1;
                while (egNorm/divisor > gradientNormBound)
                    %disp('Inflated divisor!');
                    divisor = divisor*10;
                end
                x_t = x_t - eta*effGrad/divisor;
                hu_iters = hu_iters+1;
                
                % Check to see if L iterations have passed (triggers hessian update)
                if (mod(hu_iters,L) == 0)
                    % Increment the number of hessian correction pairs
                    r = r + 1;
                    hu_iters = 0;
                    % Finish computing u_r
                    u_r = u_r/L;
                    
                    % Use HVP to compute hessian updates. Begin by sampling a minibatch
                    sidx = randsample(1:N, bH);
                    % Compute s_r update
                    s_r = u_r - u_r_prev;
                    % Compute y_r estimate using HVP
                    [~, grad_up] = ff((u_r + fdHVPStepSize*s_r)', sidx);
                    [~, grad_dn] = ff((u_r - fdHVPStepSize*s_r)', sidx);
                    y_r = (grad_up'-grad_dn')/(2*fdHVPStepSize);
                    % Store latest values of s_r and y_r
                    add(s_r, y_r);
                    
                    % Resetting u_r for next evaluation
                    u_r_prev = u_r;
                    u_r = zeros(1, d);
                end    
            end
            
            Phi = ind((j*l*bs+1):(j*l*bs+l*bs));         % Line 13
            xBar = sum(x_t_pos, 1)/S_l;         % Line 12
            % New theta_i: Line 17; start by setting theta_m_binding
            theta_m1_binding = theta_m_binding;
            theta_m1_binding(Phi) = max(theta_m1_binding)+1;
            min_Tidx = min(theta_m1_binding);
            max_Tidx = max(theta_m1_binding);
            % Reallocate the theta matrix
            theta_m1 = zeros((max_Tidx-min_Tidx+1), d);
            for curr_Tidx = min_Tidx:(max_Tidx-1)
                theta_m1(curr_Tidx-min_Tidx+1,:) = theta_m(curr_Tidx,:);
            end
            theta_m1(max_Tidx-min_Tidx+1,:) = xBar;
            % Need to update full gradient estimate: Line 18
            [~, new_phi_grad] = ff(xBar', Phi);
            % Compute update components from previous theta_m's
            old_bindings = theta_m_binding(Phi);
            old_grad = 0;
            for curr_theta = unique(old_bindings)
                [~, og] = ff(theta_m(curr_theta,:)', Phi(old_bindings==curr_theta));
                old_grad = old_grad + sum(old_bindings==curr_theta)*og;
            end
            aBar_m = aBar_m + (new_phi_grad*length(Phi) - old_grad)/N;
            % Now adjust theta_m_binding indices
            theta_m_binding = theta_m1_binding-min(theta_m1_binding)+1;
            theta_m = theta_m1;
            
            if (norm(aBar_m)/max([1 norm(xBar)]) < epsilon)
                disp("Convergence!");
                return;
            end
        end
    end

    function [r_2] = twoLoopRecursion(v_t)
        alphas = zeros(1, currDepth);

        % Begin by cloning the input gradient
        q = v_t;

        % The first loop (starts from the latest entry and goes to the earliest)
        for curr_idx = 1:currDepth
            % Compute and store alpha_i = rho_u*s_i*q
            alpha = rho(curr_idx)*dot(s(curr_idx,:), q);
            alphas(curr_idx) = alpha;
            % Update q: q = q - alpha_i*y_i
            q = q - alpha*y(curr_idx,:);
        end
        % Start computing R; Begin by computing gamma_k = s_k*y_k/(y_k*y_k)
        gamma = dot(s(currDepth,:), y(currDepth,:))/dot(y(currDepth,:), y(currDepth,:));
        % r = gamma_k*q/(1 + delta*gamma_k); the denominator includes the 
        % pseudo-hessian regularization parameter delta NOTE: There is no 
        % need to multiply by I here, as it will anyway be a dot product 
        r_2 = q*gamma/(1.0 + delta*gamma);

        % Second loop (goes in reverse, starting from the earliest entry)
        for curr_idx = flip(1:currDepth)
            % beta = rho_i*y_i*r
            beta = rho(curr_idx)*dot(y(curr_idx,:), r_2);
            % r = r + s_i*(alpha_i-beta)
            r_2 = r_2 + s(curr_idx,:)*(alphas(curr_idx)-beta);
        end
    end

    function add(inS, inY)
        % Compute rho in the lbfgs two-loop method: rho_j = 1/s_j^T*y_j
        rho = cycleDownVector(rho, 1.0/dot(inS, inY));
        % Add values to structure
        s = cycleDownMatrix(s, inS);
        y = cycleDownMatrix(y, inY);
        currDepth = min(currDepth+1, M);
    end

    % Cycle rows in matrix downwards and add a new row on top
    function [input] = cycleDownMatrix(input, newRow)
        for i = flip(2:size(input, 1))
            input(i,:) = input(i-1,:);
        end
        input(1,:) = newRow;
    end

    function [input] = cycleDownVector(input, newValue)
        for i = flip(2:length(input))
            input(i) = input(i-1);
        end
        input(1) = newValue;
    end
end