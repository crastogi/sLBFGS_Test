function [f, df, ddf] = svm(x, varargin)
    global data;
    global classes;
    lambda = .001;
    nDataPoints = length(classes);
    
    % -- function value ---------------------------------------------------
    % vector-compute alpha value
    alpha = sum(data.*x, 2);
    % Vector-compute temp value
    temp = max(0, 1-classes.*alpha);
    f = sum(temp.^2)/(2*nDataPoints) + lambda*(x*x')/2;
    df = NaN;
    ddf = NaN;
end

