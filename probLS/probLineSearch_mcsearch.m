% Incorporates mcsearch with the superior gradient variance estimator and
% the adaptive noise concept
function [outs, alpha0_out, y_tt, dy_tt, x_tt, var_f_tt, var_df_tt] = ...
        probLineSearch_mcsearch(func, x0, f0, df0, search_direction, ...
        alpha0, verbosity, outs, paras, var_f0, var_df0, grad_matrix, mu_k, w_k)
% probLineSearch.m -- A probabilistic line search algorithm for nonlinear
% optimization problems with noisy gradients. 
%
% == INPUTS ===============================================================
% [f,f', var_f, var_df] = func(x) -- function handle 
% input: 
%     x -- column vectors (positions) (Dx1) 
% output: 
%     f -- scalar function values
%     df -- vector gradients (Dx1)
%     var_f -- estimated noise for function values (scalar)
%     var_df -- estimated noise for gradients (Dx1)
% x0 -- current position of optimizer (Dx1)
% f0 -- function value at x0 (scalar, previous output y_tt)
% df0 -- gradient at x0 ((Dx1), previous output dy_tt)
% search_direction -- (- df(x0) does not need to be normalized)
% alpha0: initial step size (scalar, previous output alpha0_out)
% var_f0 -- variance of function values at x0. (scalar, previous output var_f_tt)
% var_df0 -- variance of gradients at x0. ((Dx1), previous output var_df_tt)
% verbosity -- level of stdout output. 
%         0 -- no output
%         1 -- steps, function values, state printed to stdout
%         2 -- plots with only function values
%         3 -- plots including Gaussian process and Wolfe condition beliefs.
% paras -- possible parameters that func needs.
% outs -- struct with collected statistics 
%
% == OUTPUTS ==============================================================
% outs -- struct including counters and statistics
% alpha0_out -- accepted stepsize * 1.3 (initial step size for next step)
% x_tt -- accepted position
% y_tt -- functin value at x_tt
% dy_tt -- gradient at x_tt
% var_f_tt -- variance of function values at x_tt
% var_df_tt -- variance of gradients values at x_tt
%
% Copyright (c) 2015 (post NIPS 2015 release 4.0), Maren Mahsereci, Philipp Hennig 
% mmahsereci@tue.mpg.de, phennig@tue.mpg.de
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the <organization> nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Declare global variables
global batchsize data;
% These global variabels are for mcstep/mcsearch
global isBracketed alphaMin alphaMax alphaU fU gU alphaL fL gL;
% -- setup fixed parameters -----------------------------------------------
if ~isfield(outs, 'counter')
    outs.counter = 1;
end
if ~isfield(outs, 'alpha_stats')
    outs.alpha_stats = alpha0; % running average over accepted step sizes
end
if isempty(verbosity)
    verbosity = 0;
end
% -- set up GP ------------------------------------------------------------
% create variables with shared scope. Ugly, but necessary because
% matlab does not discover shared variables if they are first created in a
% nested function.
d2m = []; d3m = []; V = []; Vd = []; dVd = [];
m0  = []; dm0 = []; V0= [];Vd0 = []; dVd0= [];
V0f = []; Vd0f= []; V0df=[]; Vd0df = [];

% kernel:
offset  = 10; % off-set, for numerical stability.
k  = @(a,b) (3 \ min(a+offset,b+offset).^3 + 0.5 * abs(a-b) .* min(a+offset,b+offset).^2);
kd = @(a,b) ((a<b) .* ((a+offset).^2/2) + (a>=b) .* ((a+offset)*(b+offset) - 0.5 .* (b+offset).^2));
dk = @(a,b) ((a>b) .* ((b+offset).^2/2) + (a<=b) .* ((a+offset)*(b+offset) - 0.5 .* (a+offset).^2));
dkd= @(a,b) (min(a+offset,b+offset));

% further derivatives
ddk = @(a,b) (a<=b) .* (b-a);
ddkd= @(a,b) (a<=b);
dddk= @(a,b) -(a<=b);

% -- helper functions -----------------------------------------------------
GaussCDF = @(z) 0.5 * (1 + erf(z/sqrt(2)));

% LINE SEARCH PARAMETERS
% constants for Wolfe conditions (must be chosen 0 < c1 < c2 < 1)
c1 = 0.0001;   
c2 = 0.99; 
WolfeThreshold = 0.3; % Condition to trigger acceptance of wolfe criteria
maxFuncEvals = 7; % maximum #function evaluations in one line search (+1)
smallStepSize = 0.01; % size of non-dimensional step to take when failure occurs
lastStep = false;
improvedVariance = false;
stepAlphaMin = 1e-20;
stepAlphaMax = 1e20;
maxMCSearchIters = 10;
uTol = 1e-10;
if ~isempty(mu_k) && ~isempty(w_k)
    % Use SVRG for gradient moves
    useSVRG = true;
else
    useSVRG = false;
end

% BEGIN LINE SEARCH
tt  = 1; % initial step size in scaled space
% Compute scaling factor beta
beta = abs(search_direction'*df0); % scale f and df according to 1/(beta*alpha0)
% Compute noise for starting point (T=0)?
sigmaf  = sqrt(var_f0)/(alpha0*beta); 
if improvedVariance
    % Compute improved variance
    g0dot   = zeros(1, size(grad_matrix, 1));
    for currIdx = 1:length(g0dot)
        g0dot(currIdx) = grad_matrix(currIdx,:)*search_direction;
    end
    sigmadf = sqrt(1/(batchsize-1)*(1/batchsize*sum(g0dot.^2)-(sum(g0dot)/batchsize)^2))/beta;
else
    sigmadf = sqrt((search_direction.^2)'*var_df0)/beta;
end

% -- initiate data storage ------------------------------------------------
% VERIFY THAT THESE NEED TO EXIST (todo)
T            = 0;
Y            = 0; 
dY           = df0;
dY_projected = (df0'*search_direction)/beta;  
Sigmaf       = var_f0;
Sigmadf      = var_df0;
N            = 1;

% Update GP at starting point
updateGP();

% Loop until budget is consumed 
while true
    % Begin by evaluating function and updating GP at new tt value
    evaluate_function();
    
    % Test for NaN or Inf and reset
    if (isnan(y) || isinf(y) || any(isnan(dy)+isinf(dy)))
        if verbosity >0
            disp('NaN Error!');
        end
        % Rescale search direction to unit vector and reset evaluations
        search_direction = search_direction/norm(search_direction);
        tt = 1;
        beta = abs(search_direction'*df0);
        sigmaf = sqrt(var_f0)/(alpha0*beta);
        g0dot   = search_direction;
        for currIdx = 1:length(grad_matrix)
            g0dot(currIdx) = grad_matrix(currIdx,:)*search_direction;
        end
        sigmadf = sqrt(1/(batchsize-1)*(1/batchsize*sum(g0dot.^2)-(sum(g0dot)/batchsize)^2))/beta;
        T = 0;
        Y = 0;
        dY = df0;
        dY_projected = (df0'*search_direction)/beta;  
        Sigmaf       = var_f0;
        Sigmadf      = var_df0;
        N            = 1;
        updateGP();
        continue;
    end
    
    % No numerical instability; update GP
    updateGP();
    
    % Next, verify that the GP gradient at the origin is negative. If it is
    % not, it means the GP believes the objective function is increasing,
    % and the gradient direction is NOT a descent direction. Take a really
    % small step and terminate.
     if d1m(0) > 0
        if verbosity > 0
            disp(['GP detects that the current direction is not a descent direction; taking a small step and terminating. Number of line search iterations: ' num2str(N) '; T string: ' num2str(T')]);
            %gp_wolfe_diag();
        end
        % Evaluate function at small step size
        %tt = smallStepSize;
        tt = -1;
        evaluate_function();
        make_outs(y, dy, var_f, var_df);
        return;
    end
    
    % The following can have two alternatives: Accept point only if it
    % satisfies WC on latest point, or accept ANY point that satisfies WC
    % the most. This seems most permissive...
    % Next, confirm if the new point satisfies the wolfe conditions
    if lastStep
        if probWolfe(tt) > WolfeThreshold % are we done?
            if verbosity > 0
                disp(['Found acceptable point. Number of line search iterations: ' num2str(N) '; T string: ' num2str(T')]); 
            end
            make_outs(y, dy, var_f, var_df);
            return;   
        end
    else
        % Loop over all evaluated positions, except the original point; use
        % Wolfe probabilities and GP means to accept points
        wolfe_values= zeros(1, length(T));      % Ensure idx=1 has WP=0
        gp_means    = wolfe_values;
        % Start by evlauating GP and wolfe probability
        for currT = 2:length(T)
            wolfe_values(currT) = probWolfe(T(currT));
            gp_means(currT)     = m(T(currT));  
        end
        % See what points match wolfe criteria
        wolfe_idx = wolfe_values>WolfeThreshold;
        if sum(wolfe_idx)>0
            % Points exist that match wolfe criteria. Find optimal one
            if verbosity > 0
                disp(['Found acceptable point. Number of line search iterations: ' num2str(N) '; T string: ' num2str(T')]); 
            end
            T_temp  = T(wolfe_idx);
            gp_means= gp_means(wolfe_idx);
            [~, min_gp_idx] = min(gp_means);
            tt = T_temp(min_gp_idx);
            
            % Find the corresponding gradient and variances of optimal pt
            dy = dY(:, T == tt); y = Y(T == tt); var_f = Sigmaf(T == tt); var_df = Sigmadf(:, T==tt);    
            make_outs(y, dy, var_f, var_df);
            return; 
        end
    end
    
    % Compute alternative acceptance criteria/see if old points satsify WC?
    % We will ignore this for now...
	% -- check this point as well for acceptance --------------------------
%     if abs(dmin) < 1e-5 && Vd(tmin) < 1e-4 % nearly deterministic
%         tt = tmin; dy = dY(:, minj); y = Y(minj); var_f = Sigmaf(minj); var_df = Sigmadf(:, minj);        
%         disp('found a point with almost zero gradient. Stopping, although Wolfe conditions not guaranteed.')
%         make_outs(y, dy, var_f, var_df);
%         return;
%     end
    
    % Check to see if the number of line search steps has been exceeded
    if N >= maxFuncEvals+1
        break;
    end
    
    % None of the currently evaluated points are acceptable. Compute new
    % test point using mcsearch.
    % (todo) Need to figure out the best way to initialize alpha for
    % mcsearch.
    try 
        tt = mcsearch(max(T));
        % Check to see if tt has already been tested (occurs when an
        % existing point satisfies the non-probabilistic wolfe conditions)
        if sum(ismember(T, tt))~=0
            %disp('tt has the same value as existing T value being tested!');
            % Check to see if tt is the same as the rightmost point tested
            T_sorted = sort(T);
            T_idx = find(T_sorted==tt);
            if (T_idx==length(T))
                % tt IS the rightmost point tested
                tt = max(T)*2;
            else 
                % Set new tt to be the average between it and the point
                % right of it
                tt = (T_sorted(T_idx)+T_sorted(T_idx+1))/2;
            end
        end
    catch mcsearch_error
        % An error has taken place in the line search. Exit making a small
        % step (todo) or unit step in the original direction
        disp('MCSEARCH FAILURE!');
        disp(mcsearch_error.message);
        disp(['T string: ' num2str(T')]);
        % Debugging test; begin by visualizing GP
        gp_wolfe_diag();
        
        % Evaluate function at small step size
        tt = smallStepSize;
        evaluate_function();
        make_outs(y, dy, var_f, var_df);
        return;
    end
end

% Algorithm reached limit without finding acceptable point; Evaluate a
% final time and return the "best" point (one with lowest function value)
% (todo): make small step or unit step?
gp_means    = inf(1, length(T));    % Ensures idx=1 has mean=inf
% Start by evlauating GP and wolfe probability
for currT = 2:length(T)
    gp_means(currT)     = m(T(currT));  
end
[~, min_gp_idx] = min(gp_means);
tt = T(min_gp_idx);
% Find the corresponding gradient and variances of optimal pt
dy = dY(:, T == tt); y = Y(T == tt); var_f = Sigmaf(T == tt); var_df = Sigmadf(:, T==tt);    
make_outs(y, dy, var_f, var_df);

% *************************************************************************
function gp_wolfe_diag()
	nSteps = range(T)*100;
    x_axis = linspace(min(T), max(T), range(T)*100);
    f_outs = zeros(1,nSteps);
    w1_vals= f_outs;
    g_outs = f_outs;
    for curr_idx = 1:nSteps 
        f_outs(curr_idx) = m(x_axis(curr_idx));
        w1_vals(curr_idx) = m(0)+x_axis(curr_idx)*c1*d1m(0);
        g_outs(curr_idx) = d1m(x_axis(curr_idx));
    end
    figure;
    yyaxis left; plot(x_axis, f_outs); hold on; plot(x_axis, w1_vals); ylabel('GP Function Value');
    yyaxis right; plot(x_axis, g_outs); ylabel('GP Gradient Value');
    hold on; line([min(x_axis) max(x_axis)], [c2*d1m(0) c2*d1m(0)]); line([min(x_axis) max(x_axis)], [-c2*d1m(0) -c2*d1m(0)]);
end

% *************************************************************************
function [alphaT] = mcsearch(alphaT)
    % Check to see if guess alphaT is within bounds
    if (alphaT<stepAlphaMin || alphaT>stepAlphaMax)
        error('Line search failure: initial step alpha guess out of bounds!');
    end
    
%    disp('Entering mcsearch.');
    
    % Initialize variables: These are trivial inputs for the GP
    %x0MC = 0;
    f0MC = m(0);
    g0MC = d1m(0);
    
    % These parameters are shared with mcstep
    alphaL			= 0;
    fL				= f0MC;
    gL				= g0MC;
    alphaU			= 0;
    fU				= f0MC;
    gU				= g0MC;
    currIntWidth	= stepAlphaMax-stepAlphaMin;
    prevIntWidth	= 2*currIntWidth;
    useModifiedFunc	= true;
    isBracketed		= false;
    nLSEvals		= 0;

    % Create first interval, and ensure that alphaT lies within the range
    % of acceptable alpha values
    alphaT			= min([max([alphaT stepAlphaMin]) stepAlphaMax]);
    alphaMin		= alphaL;
    alphaMax		= alphaT + 4*(alphaT-alphaL);

    while nLSEvals <= maxMCSearchIters
        
%        disp('Executing one line search step.');
        
        % Evaluate function
        fT	= m(alphaT);
        gT	= d1m(alphaT);
        nLSEvals = nLSEvals + 1;

        % Check for convergence
        fTTest = f0MC+alphaT*c1*g0MC;
        if (alphaT==stepAlphaMax && (fT-fTTest<=0 && gT-c1*g0MC<0))
            % \psi(\alpha_t) \leq 0 and \psi^'(\alpha_t) < 0
            error('Line search failure: search terminated, step alpha at maximum!');
        elseif (alphaT==stepAlphaMin && (fT-fTTest>0 || gT-c1*g0MC>=0))
            % \psi(\alpha_t) > 0 and \psi^'(\alpha_t) \geq 0
            error('Line search failure: search terminated, optimal step alpha is lower than minimum!');
        end 
        % Converged under strong Wolfe conditions
        if (fT-fTTest<=0 && abs(gT)<=c2*abs(g0MC))
            return;
        end

        % Check if modified update should be used: 
        if (useModifiedFunc && fT-fTTest<=0 && gT>=0)
            % \psi(\alpha_t) \leq 0 and \phi^'(\alpha_t)>0
            useModifiedFunc = false;
        end
        % Generate a new safeguarded alphaT and interval I
        if useModifiedFunc
            alphaT = mcstep(alphaL, fL-alphaL*c1*g0MC, gL-c1*g0MC, alphaT, fT-alphaT*c1*g0MC, gT-c1*g0MC, alphaU, fU-alphaU*c1*g0MC, gU-c1*g0MC);
        else
            alphaT = mcstep(alphaL, fL, gL, alphaT, fT, gT, alphaU, fU, gU);
        end
        
        % Ensure that the interval I decreases sufficiently using a bisection prediction  
        if (isBracketed && (abs(alphaU-alphaL)>=.66*prevIntWidth))
            alphaT = alphaL + (alphaU-alphaL)/2;
        end
        alphaT = min([max([alphaT stepAlphaMin]) stepAlphaMax]);
        % Define bounds of new interval
        if isBracketed
            alphaMin = min([alphaL alphaU]);
            alphaMax = max([alphaL alphaU]);
        else
            alphaMin = alphaL;
            alphaMax = alphaT+4*(alphaT-alphaL);
        end
        prevIntWidth = currIntWidth;
        currIntWidth = alphaMax-alphaMin;
        % Check to see that the interval is larger than uTol (machine precision)
        if (isBracketed && currIntWidth<=uTol*alphaMax)
            error('Line search failure: The relative width of the interval of uncertainty is less than uTol!');
        end
        % If unusual termination is about to occur, let alphaT be the lowest point found so far
        if ( (isBracketed && (alphaT<=alphaMin || alphaT>=alphaMax)) || nLSEvals>=maxMCSearchIters-1) 
            if verbosity > 1
                disp("UNUSUAL TERMINATION!");
            end
            alphaT = max(T)*2;
            %gp_wolfe_diag();
            return;
        end
    end
    error('Line search failure: number of line search function calls exceeded maximum!');
end

% *************************************************************************
function evaluate_function()

    outs.counter = outs.counter + 1;
    if useSVRG
        sidx = randsample(1:(size(data,1)), batchsize);
        [y, dy, var_f, var_df,~,~,grad_matrixT] = func(x0 + tt*alpha0*search_direction, sidx);
        [~,grad_f_wk,~,~,~,~] = func(w_k', sidx);
        dy = (dy + mu_k') - grad_f_wk;
    else
        % y: function value at tt
        [y, dy, var_f, var_df,~,~,grad_matrixT] = func(x0 + tt*alpha0*search_direction);
    end
    
    if isinf(y) || isnan(y)
        % this does not happen often, but still needs a fix
        % e.g. if function return inf or nan (probably too large step), 
        % do function value clipping relative to the intial value, 
        % e.g. y = 1e3*f0. 
        error('function values is inf or nan.')
    end
    
    % -- scale f and df ---------------------------------------------------
    y            = (y - f0)/(alpha0*beta);        % subtract offset    
    dy_projected = (dy'*search_direction)/beta;   % projected gradient 
    
    % -- store ------------------------------------------------------------
    T            = [T; tt]; 
    Y            = [Y; y]; 
    dY           = [dY, dy];
    dY_projected = [dY_projected; dy_projected]; 
    Sigmaf       = [Sigmaf; var_f];
    Sigmadf      = [Sigmadf, var_df];
    N            = N + 1;
    
    % Take max of observed sigmaf, sigmadf
    sigmaf_test  = sqrt(var_f)/(alpha0*beta);
    if improvedVariance
        g0dot2   = zeros(1, size(grad_matrixT, 1));
        for currIdx2 = 1:length(g0dot2)
            g0dot2(currIdx2) = grad_matrixT(currIdx2,:)*search_direction;
        end
        sigmadf_test = sqrt(1/(batchsize-1)*(1/batchsize*sum(g0dot2.^2)-(sum(g0dot2)/batchsize)^2))/beta;
    else
        sigmadf_test = sqrt((search_direction.^2)'*var_df)/beta;
    end
    
    if sigmaf_test > sigmaf
        sigmaf = sigmaf_test;
    end
    if sigmadf_test > sigmadf
        sigmadf = sigmadf_test;
    end
end

% -- helper functions -----------------------------------------------------
function updateGP() % using multiscope variables to construct GP

    % build Gram matrix
    kTT   = zeros(N); 
    kdTT  = zeros(N); 
    dkdTT = zeros(N);
    for i = 1:N
        for j = 1:N
            kTT(i,j)   = k(T(i),  T(j));
            kdTT(i,j)  = kd(T(i), T(j));
            dkdTT(i,j) = dkd(T(i),T(j));
        end
    end
    
    % build noise matrix
    Sig = sigmaf^2 * ones(2*N, 1); Sig(N+1:end) = sigmadf^2;
    
    % build Gram matrix
    G = diag(Sig) + [kTT, kdTT; kdTT', dkdTT];
    A = G \ [Y; dY_projected];

    % posterior mean function and all its derivatives
    m   = @(t) [k(t, T')   ,  kd(t,  T')] * A;
    d1m = @(t) [dk(t, T')  , dkd(t,  T')] * A;
    d2m = @(t) [ddk(t, T') ,ddkd(t,  T')] * A;
    d3m = @(t) [dddk(t, T'),zeros(1, N)]  * A;
    
    % posterior marginal covariance between function and first derivative
    V   = @(t) k(t,t)   - ([k(t, T') ,  kd(t, T')] * (G \ [k(t, T') , kd(t, T')]'));
    Vd  = @(t) kd(t,t)  - ([k(t, T') ,  kd(t, T')] * (G \ [dk(t, T'),dkd(t, T')]'));
    dVd = @(t) dkd(t,t) - ([dk(t, T'), dkd(t, T')] * (G \ [dk(t, T'),dkd(t, T')]'));
       
    % belief at starting point, used for Wolfe conditions
    m0   = m(0);
    dm0  = d1m(0);    
    V0   = V(0);
    Vd0  = Vd(0);
    dVd0 = dVd(0);
    
    % covariance terms with function (derivative) values at origin
    V0f   = @(t) k(0,t)  - ([k(0, T') ,  kd(0, T')] * (G \ [k(t, T') , kd(t, T')]'));
    Vd0f  = @(t) dk(0,t) - ([dk(0, T'), dkd(0, T')] * (G \ [k(t, T') , kd(t, T')]'));
    V0df  = @(t) kd(0,t) - ([k(0, T'),   kd(0, T')] * (G \ [dk(t, T'),dkd(t, T')]'));
    Vd0df = @(t) dkd(0,t)- ([dk(0, T'), dkd(0, T')] * (G \ [dk(t, T'),dkd(t, T')]'));
end

function [p,p12] = probWolfe(t) % probability for Wolfe conditions to be fulfilled

    % marginal for Armijo condition
    ma  = m0 - m(t) + c1 * t * dm0;
    Vaa = V0 + (c1 * t).^2 * dVd0 + V(t) + 2 * (c1 * t * (Vd0 - Vd0f(t)) - V0f(t));
    
    % marginal for curvature condition
    mb  = d1m(t) - c2 * dm0;
    Vbb = c2^2 * dVd0 - 2 * c2 * Vd0df(t) + dVd(t);
    
    % covariance between conditions
    Vab = -c2 * (Vd0 + c1 * t * dVd0) + V0df(t) + c2 * Vd0f(t) + c1 * t * Vd0df(t) - Vd(t);
                                     
    if (Vaa < 1e-9) && (Vbb < 1e-9) % deterministic evaluations
        p = (ma >= 0) .* (mb >= 0);
        return
    end
    
    % joint probability 
    rho = Vab / sqrt(Vaa * Vbb);
    if Vaa <= 0 || Vbb <= 0
        p   = 0; 
        p12 = [0,0,0]; 
        return;
    end
    upper = (2 * c2 * (abs(dm0)+2*sqrt(dVd0))-mb)./sqrt(Vbb);
    p = bvn(-ma / sqrt(Vaa), inf, -mb / sqrt(Vbb), upper, rho);
    
    if nargout > 1
        % individual marginal probabilities for each condition 
        % (for debugging)
        p12 = [1 - GaussCDF(-ma/sqrt(Vaa)), ....
            GaussCDF(upper)-GaussCDF(-mb/sqrt(Vbb)),...
            Vab / sqrt(Vaa * Vbb)];
    end
end

function make_outs(y, dy, var_f, var_df)
    
    x_tt      = x0 + tt*alpha0*search_direction; % accepted position
    y_tt      = y*(alpha0*beta) + f0;            % function value at accepted position
    dy_tt     = dy;                              % gradient at accepted position
    var_f_tt  = var_f;                           % variance of function value at accepted position
    var_df_tt = var_df;                          % variance of gradients at accepted position
        
    % Return output step size
    outs.step_size = tt*alpha0;
    
    % (todo) get rid of all these alpha resets
    % set new set size
    % next initial step size is 1.3 times larger than last accepted step size
    alpha0_out = tt*alpha0 * 1.3;

    % running average for reset in case the step size becomes very small
    % this is a saveguard
    gamma = 0.95;
    outs.alpha_stats = gamma*outs.alpha_stats + (1-gamma)*tt*alpha0;

    % reset NEXT initial step size to average step size if accepted step
    % size is 100 times smaller or larger than average step size
    if (alpha0_out > 1e2*outs.alpha_stats)||(alpha0_out < 1e-2*outs.alpha_stats)
        if verbosity > 0
            disp 'making a very small step, resetting alpha0'; 
        end
        alpha0_out = outs.alpha_stats; % reset step size
    end
end
end