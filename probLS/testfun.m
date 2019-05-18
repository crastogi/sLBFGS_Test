function [f, df, ddf] = testfun(x, varargin)
    global data;
    nDataPoints = size(x, 1);
    f = zeros(nDataPoints,1);
    for idx = 1:nDataPoints
        f(idx) = sum((data(:,1)-x(idx,1)).^2)+sum((data(:,2)-x(idx,2)).^2);
    end
    f = f/size(data,1);
    df = NaN;
    ddf = NaN;
end

