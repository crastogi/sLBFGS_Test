function [f, df, ddf] = mccormick2(x, varargin)
% McCormick function
% usage [f, df, ddf] = mccormick2(x, varargin)
% 
% N -- # datapoints (points to evaluate)
% x   in R^[Nx2]   
% f   in R^[Nx1]   function value
% df  in R^[Nx2]   gradient 
% ddf in R^[Nx2x2] hessian
%
% (C) 2013, Maren Mahsereci (mmahsereci@tue.mpg.de)

    % -- helper parameters ------------------------------------------------ 
    a = sin(x(:,1) + x(:,2));
    b = cos(x(:,1) + x(:,2));
    
    % -- function value ---------------------------------------------------
    f =   a + (x(:,1) - x(:,2)).^2 - 1.5.*x(:,1) + 2.5.* x(:,2) + 1;
    
    % -- gradient ---------------------------------------------------------
    if nargout > 1
        df = [b + 2*(x(:,1) - x(:,2)) - 1.5, ...
              b - 2*(x(:,1) - x(:,2)) + 2.5];
    end
    
    % -- hessian ----------------------------------------------------------
    if nargout > 2
        ddf = zeros(size(x,1),2,2);
        ddf(:,1,1) = -a +2;
        ddf(:,1,2) = -a -2;
        ddf(:,2,1) = ddf(:,1,2);
        ddf(:,2,2) = ddf(:,1,1);
    end
end

