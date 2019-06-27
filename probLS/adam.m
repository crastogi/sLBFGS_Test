% Reset data consumption
global nDataPoints
nDataPoints = 0;

alpha = 0.005;                          % Stepsize
b1 = .9;                                % Beta 1 parameter
b2 = .999;                              % Beta 2 parameter
e = 1E-8;                               % Epsilon
bs = 20;                                % Batchsize
M = 1000;                               % M: max epochs
N = size(data,1);                       % Number of data points
d = size(data,2);                       % Number of dimensions
epsilon = 1E-5;                         % Convergence criteria
nIters = floor(N/bs);                   % Number of iterations per epoch

x0 = zeros(1,d);                        % Start position
xt = x0;
wk = x0;
mt = 0;
vt = 0;
totIters = 0;
sums = 0;
for m = 1:M                             % Loop over epochs
    [~, muk] = func(wk', 1:N);          % SVRG update
    disp([num2str(log10(norm(muk))) ' ' num2str(log10(norm(xt-x_min{1}')))]);
    if (norm(muk)/max([1 norm(wk)]) < epsilon)
        disp("Convergence!");
        return;
    end
    sums = xt;
    for t = 1:nIters                    % Loop over iterations
        totIters = totIters+1;
        idx = randsample(1:N, bs);
        [~, gt] = func(xt', idx);
        [~, gk] = func(wk', idx);
        gt = gt - gk + muk;
        mt = b1*mt + (1-b1)*gt';
        vt = b2*vt + (1-b2)*(gt.^2)';
        mHat = mt/(1-b1^totIters);
        vHat = vt/(1-b2^totIters);
        xt = xt - alpha*mHat/(sqrt(sum(vHat))+e);
        sums = sums + xt;
    end
    sums = sums - xt;
    wk = sums/nIters;
end

%     disp([num2str(log10(norm(aBar_m))) ' ' num2str(log10(norm(xBar'-x_min{1}))) ' '...
%num2str(sum(theta_m_binding==1)*(sum(theta_m(1,:))==0))]);
%     if (norm(aBar_m)/max([1 norm(xBar)]) < epsilon)
%         disp("Convergence!");
%         return;
%     end