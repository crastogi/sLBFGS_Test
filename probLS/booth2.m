function [f, df, ddf] = booth2(x, varargin)
% booth's function
% usage [f, df, ddf] = booth2(x, varargin)
% 
% N -- # datapoints (points to evaluate)
% x   in R^[Nx2]   
% f   in R^[Nx1]   function value
% df  in R^[Nx2]   gradient 
% ddf in R^[Nx2x2] hessian
%
% (C) 2013, Maren Mahsereci (mmahsereci@tue.mpg.de)


    % -- helper parameters ------------------------------------------------ 
    a = x(:,1) + 2*x(:,2) -7;
    b = x(:,2) + 2*x(:,1) -5;
    
    % -- function value ---------------------------------------------------
    f =   a.^2 + b.^2;
    
    % -- gradient ---------------------------------------------------------
    if nargout > 1
        df = [2*a + 4*b, 4*a + 2*b];
    end
    
    % -- hessian ----------------------------------------------------------
    if nargout > 2
        ddf = zeros(size(x,1),2,2);
        ddf(:,1,1) = 10;
        ddf(:,1,2) = 8;
        ddf(:,2,1) = ddf(:,1,2);
        ddf(:,2,2) = 10;
    end
end

