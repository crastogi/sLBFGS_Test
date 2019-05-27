% Minimal example for running the sgd+probabilistic line search on noisy
% functions. 
%
% The line search has access to exact values of sigmaf and sigmadf
%
% To get an intuition you can change the following settings
%   testfunction -- choose among 3 test functions (1, 2, 3)
%   T -- # function evaluations
%   sigmaf -- absolute noise of function values
%   sigmadf -- absolute noise of gradient evaluations 
%   verbosity -- 0: silence, 1: speak, 2: plot simple, 3: plot full (slow 
%     but good for intuition)
% and run several times. 
% 
% In this synthetic setting the noise sigmaf, sigmadf are absolute. 
% This means that the signal-to-noise ratio gets worse with dropping 
% function and gradient value (mostly closer to the minimum). 
% This also means diffusion towards the end. This is clearly visible 
% in the optimization path. Noise in real world applicatons might depend 
% on the location x, thus they need to be estimated locally.
%
% (C) 2015, Maren Mahsereci (mmahsereci@tue.mpg.de)

clearvars -except totEpochs nSamples nDataPoints totDist;
global path nEpochs x_min data batchsize stochIters;
global vf vdf;

% == CHANGE STUFF BELOW HERE ==============================================
verbosity    = 1; % 0: silence, 1: speak, 2: plot simple, 3: plot full (slow)
ff           = @noisyFunction;
maxEpochs    = 1000;
stochIters   = 5;
testfunction = 5;  % choose among 3 test functions (1, 2, 3)
probLSFunc   = @probLineSearch_mcsearch;
suppressPlot = 0;

% sythetic noise standard deviations for function value and gradient
sigmaf  = .01;
sigmadf = .01*ones(2, 1);
% == CHANGE STUFF ABOVE HERE ==============================================

switch testfunction
    case 1 % Booth2
        x1min = -10; x1max = 10; % x-limits for plotting
        x2min = -10; x2max = 10; % y-limits for plotting
        F        = @booth2;      % function handle
        x_min{1} = [1, 3];       % minimizer
        paras    = [];           % function does not have parameters
        x0       = [4; -5];      % startvalue
        ff(F, sigmaf, sigmadf);        % set up function 

    case 2 % McCormick2
        x1min = -1.5; x1max = 4;       % x-limits for plotting
        x2min = -3;   x2max = 4;       % y-limits for plotting
        F        = @mccormick2;        % function handle
        x_min{1} = [-0.5472, -1.5472]; % minimizer 1
        x_min{2} = [2.5944, 1.5944];   % minimizer 2       
        paras    = [];                 % function does not have parameters
        x0       = [3; -2];            % startvalue
        ff(F, sigmaf, sigmadf);        % set up function
        
    case 3 % Branin2
        x1min = -5; x1max = 10;      % x-limits for plotting
        x2min = 0; x2max = 15;       % y-limits for plotting
        F        = @branin2;         % function handle
        x_min{1} = [pi, 2.275];      % minimizer 1
        x_min{2} = [-pi,12.275];     % minimizer 2
        x_min{3} = [9.42478, 2.475]; % minimizer 3
        paras    = [];               % function does not have parameters
        x0       = [-2; 7];          % startvalue
        ff(F, sigmaf, sigmadf);        % set up function 
    
    case 4 % New Test function
        disp("Gaussian test function");
        data        = importdata("bivariate_normal.tsv");
        data        = data(1:1000,:);
        batchsize   = 10;
        clear x_min;
        F           = @testfun;
        x_min{1}    = mean(data)';          % Only minimizer
        x1min = -10; x1max = 10;        % x-limits for plotting
        x2min = -10; x2max = 10;        % y-limits for plotting
        ff          = @testfun_full;
        x0          = [-10;8];
        paras       = [];
        
    case 5 % New Branin Test function
        disp("Branin test function");
        data        = importdata("braninsample.tsv");
        data        = data(1:10000,:);
        batchsize   = 5000;
        clear x_min;
        F           = @branin_testfun;
        x0          = [-10;8];
        %x_min{1}    = pattern_search(x0', F);      % minimizer 1
        x1min = -10; x1max = 10;      % x-limits for plotting
        x2min = -10; x2max = 10;       % y-limits for plotting
        ff          = @branin_testfun_full;
        paras       = [];
    case 6 % svm
        disp("SVM Test");
        global classes;
        data        = importdata("../java/slbfgs/src/featurized_df.csv");
        data        = data.data;
        classes     = data(:,end);
        data        = data(:,1:(end-1));
        batchsize   = 200;
        clear x_min;
        F           = @svm;
        x_min{1}    = importdata("svm_min.csv")';
        x0          = flip(x_min{1})*5;
        ff          = @svm_full;
        paras       = [];
end % switch

nEpochs = 0;
xt = x0;
outs.counter        = 0;
full_dataset_iters   = 0;
alpha0 = 1;

% Storage only for plot
path                = x0;
grad_norm           = double.empty;
function_values     = ff(x0);

% Loop over epochs
while nEpochs < maxEpochs
    % Compute full gradient for variance reduction
    temp_batchsize = batchsize;
    batchsize = length(data);
    % Set local function and gradient variance
    wk = xt;
    [fFull, muk, vf, vdf, ~, ~] = ff(wk);
    % (Potentially) rescale full batch noise
    %vf = vf*min(1, norm(muk))^2*0;
    %vdf = vdf*min(1, norm(muk))^2*0;
    full_dataset_iters = full_dataset_iters + 1;
    batchsize = temp_batchsize;
    
    % Check for convergence: Exit on estimation on full dataset
    if (norm(muk)/max(1,norm(wk)) < 1E-3)
        disp("Converged!");
        break
    end
    
    % Perform stochastic iterations
    for t = 1:stochIters
        % Compute stochastic gradient estimate
        [f_xt,grad_f_xt,~,~,var_f,var_df, g_m] = ff(xt); %,g_m
        [f_wk,grad_f_wk,~,~,~,~] = ff(wk);
        df = grad_f_xt - grad_f_wk + muk;
        outs.counter = outs.counter + 2;
        
        % Update search direction
        search_direction = -df;
        
        % Line search
%        alpha0 = 1;   % fix
        [outs, alpha0, f, df, xt, var_f, var_df] = probLSFunc(ff, xt, f_xt, df, search_direction, alpha0, verbosity, outs, paras, var_f, var_df, g_m);
        
        % -- storage only for plot --------------------------------------------
        path            = [path, xt];           %#ok<AGROW>
        grad_norm       = [grad_norm, norm(muk)];
        function_values = [function_values, f]; %#ok<AGROW>
    end
    nEpochs = nEpochs + 1;
end


%% -- plot -----------------------------------------------------------------

if suppressPlot==0 
    clf;
    if (testfunction<6) 
        ax = linspace(x1min, x1max, 200);
        ay = linspace(x2min, x2max, 200);
        [X, Y] = meshgrid(ax, ay); % a x b
        XX = [X(:),Y(:)]; % a*b x 2
        Z = F(XX);
        Z = reshape(Z,size(X));
        subplot(1, 1, 1); hold on;
        if min(Z(:)) < 0
            Z = Z - min(Z(:));
        end

        gray_r = 1-gray;
        colormap(gray_r)
        imagesc([x1min, x1max],[x2min, x2max], log(Z)); axis tight; colorbar;
        plot(path(1,:), path(2,:), '-gx')
        for i=1:numel(x_min)
            plot(x_min{i}(1), x_min{i}(2), 'ro', 'markersize', 10)
        end
        figure
        dist = sqrt((path(1,:)-x_min{1}(1)).^2+(path(2,:)-x_min{1}(2) ).^2);
        plot(flip(log10(dist)), 1:length(dist));
    else
        % SVM case; just plot how convergence is occuring as a function of
        % epochs
        dist2opt = path;
        for i = 1:length(x_min{1})
            dist2opt(i,:) = dist2opt(i,:)-x_min{1}(i);
        end
        dists = sum(dist2opt.^2, 1);
        yyaxis left;
        plot(log10(sqrt(dists)));
        yyaxis right;
        plot(grad_norm);
    end
end