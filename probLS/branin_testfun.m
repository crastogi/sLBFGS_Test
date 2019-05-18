function [f, df, ddf] = testfun(x, varargin)
    global data;
    nDataPoints = size(x, 1);
    f = zeros(nDataPoints,1);
    
    a = 1;
    b = 5.1/(4*pi*pi);     % Need to estimate these
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    
    % -- function value ---------------------------------------------------
    for i = 1:nDataPoints
        A = (data(:,2)+x(i,2)) - b*(data(:,1)+x(i,1)).^2 + c*(data(:,1)+x(i,1)) - r;
        f(i) = sum(a*A.^2 + s*(1-t)*cos(data(:,1)+x(i,1)) + 10);
    end
    f = f/size(data,1);
    df = NaN;
    ddf = NaN;
end

