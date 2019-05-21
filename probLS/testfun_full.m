function [f, df, vf2, vdf2, var_f, var_df, grad_matrix] = testfun_full(x, varargin)
    global data;
    global batchsize;
    global nDataPoints;
    sampleidx = randsample(1:(size(data, 1)), batchsize);
    f = sum((data(sampleidx,1)-x(1)).^2+(data(sampleidx,2)-x(2)).^2)/batchsize;
    df = [sum(data(sampleidx,1)-x(1)); sum(data(sampleidx,2)-x(2))]*-2/batchsize;
    var_f = 1/(batchsize-1)*(1/batchsize*sum(((data(sampleidx,1)-x(1)).^2+(data(sampleidx,2)-x(2)).^2).^2)-f^2);
    var_df_x = 1/batchsize*sum((-2*(data(sampleidx,1)-x(1))).^2)-df(1)^2;
    var_df_y = 1/batchsize*sum((-2*(data(sampleidx,2)-x(2))).^2)-df(2)^2;
    var_df = 1/(batchsize-1)*[var_df_x; var_df_y];

%    global vf;
%    global vdf;
%    vf2 = vf*min(1, norm(df))^2; %/10000;
%    vdf2 = vdf*min(1, norm(df))^2; %/10000;
    vf2 = var_f;
    vdf2 = var_df;
    grad_matrix = [data(sampleidx, 1)-x(1) data(sampleidx, 2)-x(2)]*-2;

    nDataPoints = nDataPoints + batchsize;
end

