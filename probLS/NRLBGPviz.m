function NRLBGPviz(gpIdx)
    % Find the number of available GPs to scan from
    nGPs = length(regexp(fileread('/Users/chaitanya/Documents/GitWorkspaces/slbfgs/java/output/GP_Path.txt'), 'NEW GP'));
    % Make sure gpIdx >= nGPs; also, this will be the distance from LAST GP
    if (gpIdx > nGPs)
        error('gpIdx is out of range, nGPs= %f', nGPs);
    end
    % Compute number of GPs to skip
    nSkips = nGPs-gpIdx;
    fid=fopen('/Users/chaitanya/Documents/GitWorkspaces/slbfgs/java/output/GP_Path.txt');
    % First 2 lines are c1, c2 values
    c1 = str2double(fgetl(fid));
    c2 = str2double(fgetl(fid));
    % offset into correct GP line
    currSkips = 0;
    isHit = false;
    tline = fgetl(fid);
    while ischar(tline)
        if strcmp(tline, '====NEW GP====');
            currSkips = currSkips + 1;
            if (currSkips==nSkips)
                isHit = true;
                break
            end
        end
        tline = fgetl(fid);
    end
    if (~isHit) 
        error('index does not exist');
    end
    % Next line is beta value
    beta = str2double(fgetl(fid));
    tline = fgetl(fid);
    while ischar(tline)
        if strcmp(tline, '====NEW GP====')
            % Reached next GP block
            break
        end
        disp(tline);
        % Parse backwards; first see if there are any function values
        containsFunction = contains(tline, 'function');
        if containsFunction
            info = strsplit(tline, '<function_val>');
            fvals = str2double(strsplit(info{2}, ','));
            fvals = fvals(~isnan(fvals));
            fvals = (fvals-fvals(1))/beta;          % Transform
            info = strsplit(info{1}, '<function_pos>');
            fpos = str2double(strsplit(info{2}, ','));
            fpos = fpos(~isnan(fpos));
            tline = info{1};
        end
        info = strsplit(tline, '<sigmadf>');
        sigmadf = str2double(info{2});
        info = strsplit(info{1}, '<sigmaf>');
        sigmaf = str2double(info{2});
        % dYproj
        info = strsplit(info{1}, '<dYproj>');
        dY_projected = str2double(strsplit(info{2}, ','));
        dY_projected = dY_projected(~isnan(dY_projected))';
        % Y
        info = strsplit(info{1}, '<Y>');
        Y = str2double(strsplit(info{2}, ','));
        Y = Y(~isnan(Y))';
        % T
        info = strsplit(info{1}, '<T>');
        T = str2double(strsplit(info{2}, ','));
        T = T(~isnan(T))';
        % N
        info = strsplit(info{1}, '<N>');
        N = str2double(info{2});

        % Now set up GP
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

        % Plot
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
        if containsFunction
            %fill([fpos;flipud(fpos)],[fvals-sigmaf;flipud(fvals+sigmaf)],[.95 .98 .98],'linestyle','none');
            plot(fpos, fvals, '-.');
        end
        yyaxis right; plot(x_axis, g_outs); ylabel('GP Gradient Value');
        hold on; 
        line([min(x_axis) max(x_axis)], [c2*d1m(0) c2*d1m(0)], 'LineStyle', '--'); 
        line([min(x_axis) max(x_axis)], [-c2*d1m(0) -c2*d1m(0)], 'LineStyle', '--');

        tline = fgetl(fid);
    end
    fclose(fid);
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