function [f, df, ddf] = branin2(x, varargin)
% branin function
% usage [f, df, ddf] = branin2(x, varargin)
% 
% N -- # datapoints (points to evaluate)
% x   in R^[Nx2]   
% f   in R^[Nx1]   function value
% df  in R^[Nx2]   gradient 
% ddf in R^[Nx2x2] hessian
%
% (C) 2013, Maren Mahsereci (mmahsereci@tue.mpg.de)


    % -- helper parameters ------------------------------------------------ 
    a = 1;
    b = 5.1/(4*pi*pi);
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    
    A = x(:,2) - b*x(:,1).^2 + c*x(:,1) - r;

    % -- function value ---------------------------------------------------
    f =   a*A.^2 + s*(1-t)*cos(x(:,1)) + 10;
    
    % -- gradient ---------------------------------------------------------
    if nargout > 1
        df = [2*a*(-2*b*x(:,1) + c).*A - s*(1-t).*sin(x(:,1)),...
            2*a*A];
    end
    
    % -- hessian ----------------------------------------------------------
    if nargout > 2
        ddf = zeros(size(x,1),2,2);
        ddf(:,1,1) = -4*a*b*A + 2*a*(-2*b*x(:,1) + c).^2 - s*(1-t)*cos(x(:,1));
        ddf(:,1,2) = 2*a*(-2*b*x(:,1) + c);
        ddf(:,2,1) = ddf(:,1,2);
        ddf(:,2,2) = 2*a;
    end
end

