function [f, df, var_f, var_df] = noisyFunction(varargin)
% udage [f, df, var_f, var_df] = noisyFunction(varargin)
%
% wrapper function to make a 2D function noisy with independent 
% Gaussian noise in every dimension. 
% It has 2 steps:
% 1) set up (no outputs): varargin = sigmaf, sigmadf 
% 2) evaluate (with outputs): varargin = x
%
% outs: f -- noisy function value (scalar) at position x
%       df -- noisy gradient (2x1) at position x
%       var_f -- noise variance of function value (scalar) at position x
%       var_df --  noise varience of gradient (2x1) at position x
%
% (C) 2015, Maren Mahsereci (mmahsereci@tue.mpg.de)

    persistent F sigmaf sigmadf
    
    % set up persistent variables
    if nargout == 0
        p = varargin; 
        F       = p{1}; % function handle of deterministic function
        sigmaf  = p{2}; % absolute noise of function value (scalar)
        sigmadf = p{3}; % absolute noise of gradients (2x1)
       
    % evaluate function
    else
        x = varargin{1};
        [f,df] = F(x'); % noise free
        f      = f  + sigmaf .* randn();
        df     = df' + diag(sigmadf) * randn(2, 1);
        
        var_f  = sigmaf.^2;
        var_df = sigmadf.^2.*ones(2, 1);
        
    end
end