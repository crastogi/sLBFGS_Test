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

clearvars -except totEpochs nSamples nDataPoints totDist q measureMoves metrics CDLCFD CDLCFI CDLIFD CDLIFI IDLCFD IDLCFI IDLIFD IDLIFI;
global path nEpochs x_min data batchsize stochIters c1 c2 WolfeThreshold;
global vf vdf;

% == CHANGE STUFF BELOW HERE ==============================================
verbosity    = 0; % 0: silence, 1: speak, 2: plot simple, 3: plot full (slow)
ff           = @noisyFunction;
testfunction = 6;  % choose among 3 test functions (1, 2, 3)
suppressPlot = 1;

% PLS + optimizer options
useSLBFGS    = true;
maxEpochs    = 1000;
batchsize    = 20;
stochIters   = 50;
probLSFunc   = @mcsearch;   % can be null ([])
useSVRG      = true;
c1           = .1;
c2           = .9;
WolfeThreshold=.2;
variance_option=0;   % 0: standard, 1: adaptive, 2: improved

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
        data        = importdata("Data/bivariate_normal.tsv");
        data        = data(1:1000,:);
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
        data        = importdata("Data/braninsample.tsv");
        data        = data(1:10000,:);
        clear x_min;
        F           = @branin_testfun;
        x0          = [-10;8];
        x_min{1}    = [3.987363338470459; -1.360818743705750];
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
        clear x_min;
        F           = @svm;
        x_min{1}    = importdata("Data/svm_min.csv")';
        x0          = flip(x_min{1})*5;
        ff          = @svm_full;
        paras       = [];
end % switch

if useSLBFGS
    stepSize = .1;
    hessianPeriod = 10;
    memorySize = 50;
    [path, function_values, grad_norm] = sLBFGS(ff, x0, stochIters, ...
        hessianPeriod, maxEpochs, stepSize, memorySize, verbosity, ...
        probLSFunc);
    nEpochs = size(path,2);
else
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
        batchsize = size(data, 1);
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
            % Compute stochastic gradient estimate; set batch
            sidx = randsample(1:(size(data,1)), batchsize);
            [f_xt,grad_f_xt,~,~,var_f,var_df, g_m] = ff(xt, sidx);
            [f_wk,grad_f_wk,~,~,~,~] = ff(wk, sidx);
            df = grad_f_xt - grad_f_wk + muk;
            outs.counter = outs.counter + 2;

            % Update search direction
            search_direction = -df;

            % Line search with SVRG option

            if useSVRG
                [outs, alpha0, f, df, xt, var_f, var_df] = probLSFunc(ff, xt, ...
                f_xt, df, search_direction, alpha0, verbosity, outs, paras, ...
                var_f, var_df, variance_option, g_m, fFull, muk, wk);
            else
                [outs, alpha0, f, df, xt, var_f, var_df] = probLSFunc(ff, xt, ...
                f_xt, df, search_direction, alpha0, verbosity, outs, paras, ...
                var_f, var_df, variance_option, g_m, [], [], []);
            end

            % -- storage only for plot --------------------------------------------
            path            = [path, xt];           %#ok<AGROW>
            grad_norm       = [grad_norm, norm(muk)];
            function_values = [function_values, f]; %#ok<AGROW>
        end
        nEpochs = nEpochs + 1;
    end
end

%% -- plot -----------------------------------------------------------------

if suppressPlot==0
    figure;
    % Convergence information
    dist2opt = path;
    dist2pt = path*0;
    for i = 1:length(x_min{1})
        dist2opt(i,:) = dist2opt(i,:)-x_min{1}(i);
    end
    for i = 2:size(path, 2)
        dist2pt(:,i) = path(:,i)-path(:,(i-1));
    end
    dists = sum(dist2opt.^2, 1);
    ptdists = sum(dist2pt.^2, 1);
    yyaxis left;
    plot(log10(sqrt(dists)), 'DisplayName','Distance To Opt');
    hold on;
    plot(log10(sqrt(ptdists)), 'DisplayName','Epoch Pt Dist');
    ylabel('log10( Value )');
    yyaxis right;
    plot(log10(grad_norm),'DisplayName','Grad Norm');
    hold on;
    plot(log10(function_values), 'DisplayName','Function Value');
    ylabel('log10( Value )');
    legend
    if (testfunction<6)
        figure;
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
    end
end