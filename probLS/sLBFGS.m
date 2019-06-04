% Understand average distance travelled, alpha, etc. with line search and
% without line search. See where the 'distance' is being accumulated
% Play with c1/c2 
% Why does the line search 'transition' at certain points? (ls evals)

% Create automatic function profiler

function [path, function_values, grad_norm] = sLBFGS(ff, x0, ...
    epochIterations, hessianPeriod, maxEpoch, stepsize, memorySize, ...
    verbosity, ls_function)
    if (isempty(ls_function))
        useLineSearch = 0;
    else 
        useLineSearch = 1;
    end
    global batchsize data;
    b = batchsize;
    d = size(data, 2);
    N = size(data, 1);
    bH = b*10;
    M = memorySize;
    eta = stepsize;
    epsilon = 1E-5;
    m = epochIterations;
    L = hessianPeriod;
    delta = 0;
    fdHVPStepSize = 5E-2;
    gradientNormBound = 100;
    % Number of currently computed Hessian correction pairs
    r = 0;
    currDepth = 0;
    % Average of path travelled in the current and previous inverse hessian updates
    u_r = zeros(1, d);
    u_r_prev = zeros(1, d);
    % Components of two-loop update
    rho = zeros(1, M);
    s   = zeros(M, d);
    y   = zeros(M, d);

    % Initialize with a seed
    w_k = x0';
    w_k_prev = w_k;
    
    % Store path, etc.
    path = [];
    function_values = [];
    grad_norm = [];

    % Loop over epochs 
    for k = 0:maxEpoch
        % Compute full gradient for variance reduction
        batchsize = N;
        [fFull, mu_k, ~, ~, fVar_f, fVar_df] = ff(w_k');
        mu_k = mu_k';
        var_scale = (N-b)/b;
        
        path = [path, w_k'];
        function_values = [function_values, fFull];
        grad_norm = [grad_norm, norm(mu_k)]; 
        
        % Print some information about the current epoch
        if (verbosity)
            if (k==0)
                disp(['Starting Function Value: ' num2str(fFull)]);
            else
                disp(['Epoch: ' num2str(k) '; Function Value: ' num2str(fFull)]);
            end
        end

        % Check for convergence, final epoch
        if (norm(mu_k)/max([1 norm(w_k)]) < epsilon)
            disp("Convergence!");
            batchsize = b;
            return;
        end
        % This allows convergence to be tested on the FINAL epoch 
        % iteration without computation
        if (k==maxEpoch)
            break;		
        end
        x_t = w_k;					% Set x_t to current value of w_k

        % Perform m stochastic iterations before a full gradient computation takes place
        for t = 1:m
            % Compute the current stochastic gradient estimate; begin by sampling a minibatch
            batchsize = b;
            sidx = randsample(1:(size(data,1)), batchsize);
            % Next, compute the reduced variance function and gradient
            [f_xt,grad_f_xt,~,~,~,~] = ff(x_t', sidx);
            [f_wk,grad_f_wk,~,~,~,~] = ff(w_k', sidx);
            f_t = (f_xt + fFull) - f_wk;
            v_t = (grad_f_xt' + mu_k) - grad_f_wk';
            
            %disp(['New function value: ' num2str(f_xt) '; New gradient norm: ' num2str(norm(grad_f_xt))]);

            if (norm(v_t)<1E-10) 
                % Potential convergence (and scale issue, break and go to
                % outer loop to check gradient on full batch)
                break
            end
            
            % Update u_r with current position
            u_r = u_r + x_t;			
            % Compute next iteration step; condition the gradient so as
            % not to produce rapid oscilations (extreme function values)			
            if (r < 1)						
                % H_0 = I Until a one hessian correction takes place 
                effGrad = v_t;
            else
                % Compute the two-loop recursion product
                effGrad = twoLoopRecursion(v_t);
            end
            
            %Bound the effective gradient update step
            egNorm = norm(effGrad);
            if ~useLineSearch
                divisor = 1;
                while (egNorm/divisor > gradientNormBound)
                    %disp('Inflated divisor!');
                    divisor = divisor*10;
                end                
            end
            
            % Should a linesearch be used?
            if useLineSearch
                % Compute eta; different methods until one hessian 
                % correction takes place
                batchsize = b;
                if (r < 1)
                    ls_out = ls_function(ff, x_t', f_t, v_t', -effGrad',...
                        eta/norm(v_t), verbosity, [], [], fVar_f*var_scale, ...
                        fVar_df*var_scale, 0, [], fFull, mu_k', w_k');
                else
                    ls_out = ls_function(ff, x_t', f_t, v_t', -effGrad',...
                        eta, 0, [], [], fVar_f*var_scale, ...
                        fVar_df*var_scale, verbosity, [], fFull, mu_k', w_k');
                end
                x_t = x_t - ls_out.step_size*effGrad;
                %disp(['k: ' num2str(k) ' t: ' num2str(t) ' r: ' num2str(r) '  norm(effGrad): ' num2str(norm(effGrad)) ' f_t: ' num2str(f_t) ' var_f: ' num2str(var_f) ' fFull: ' num2str(fFull) ' fVar_f: ' num2str(fVar_f)]);
                disp(['k: ' num2str(k) ' t: ' num2str(t) ' r: ' num2str(r) '  norm(effGrad): ' num2str(norm(effGrad)) ' nLSEvals: ' num2str(ls_out.nLSEvals) ' stepsize: ' num2str(ls_out.step_size) ' f0: ' num2str(ls_out.f0) ' sigmaf: ' num2str(ls_out.sigmaf) ' df0: ' num2str(ls_out.df0) ' sigmadf: ' num2str(ls_out.sigmadf) ' beta: ' num2str(ls_out.beta)]);
            else
                x_t = x_t - eta*effGrad/divisor;
                %disp(['k: ' num2str(k) ' t: ' num2str(t) ' r: ' num2str(r) '  norm(effGrad): ' num2str(norm(effGrad)) ' divisor: ' num2str(divisor) ' prod: ' num2str(norm(1/norm(v_t)*effGrad))]);
            end
            
            % Check to see if L iterations have passed (triggers hessian update)
            if (mod(t,L) == 0)
                % Increment the number of hessian correction pairs
                r = r + 1;
                % Finish computing u_r
                u_r = u_r/L;

                % Use HVP to compute hessian updates. Begin by sampling a minibatch
                batchsize = bH;
                sidx = randsample(1:(size(data,1)), batchsize);
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
        % Choose either the last position vector x_t
        w_k = x_t;
    end
    batchsize = b;
    disp('sLBFGS Failure: maximum epochs exceeded without convergence!');
    return;

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