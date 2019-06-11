function [f, df, vf2, vdf2, var_f, var_df, grad_matrix, temp] = branin_testfun_full(x, sampleidx, varargin)
    a = 1;
    b = 5.1/(4*pi*pi);     % Need to estimate these
    c = 5/pi;
    r = 6;
    s = 10;
    t = 1/(8*pi);
    
    global data;
    global batchsize;
    if ~exist('sampleidx','var')
        sampleidx = randsample(1:(size(data, 1)), batchsize);
    end
    
    % -- function value ---------------------------------------------------
    A = (data(sampleidx,2)+x(2)) - b*(data(sampleidx,1)+x(1)).^2 + c*(data(sampleidx,1)+x(1)) - r;
    f = sum(a*A.^2 + s*(1-t)*cos(data(sampleidx,1)+x(1)) + 10)/batchsize;

    temp = a*A.^2 + s*(1-t)*cos(data(sampleidx,1)+x(1)) + 10;
    
    df = [sum(2*a*(c-2*b*(data(sampleidx,1)+x(1))).*A-s*(1-t)*sin(data(sampleidx,1)+x(1))); sum(2*a*A)]/batchsize;
    
    grad_matrix = [2*a*(c-2*b*(data(sampleidx,1)+x(1))).*A-s*(1-t)*sin(data(sampleidx,1)+x(1)) 2*a*A];
    
    var_f = 1/(batchsize-1)*(1/batchsize*sum((a*A.^2 + s*(1-t)*cos(data(sampleidx,1)+x(1)) + 10).^2)-f^2);
    var_df_x = 1/batchsize*sum((2*a*(c-2*b*(data(sampleidx,1)+x(1))).*A-s*(1-t)*sin(data(sampleidx,1)+x(1))).^2)-df(1)^2;
    var_df_y = 1/batchsize*sum((2*a*A).^2)-df(2)^2;
    var_df = 1/(batchsize-1)*[var_df_x; var_df_y];

    global vf;
    global vdf;
%     vf2 = vf/10000;
%     vdf2 = vdf/10000;
    vf2 = var_f; %var_f*min(1, norm(df))^3; %ceil(-log(norm(df))+1); %*(size(data,1)-batchsize)/size(data,1)*min(1, norm(df))^ceil(-log(norm(df))+1);
    vdf2 = var_df; %var_df*min(1, norm(df))^3; %ceil(-log(norm(df))+1); %(size(data,1)-batchsize)/size(data,1)*min(1, norm(df))^ceil(-log(norm(df))+1);
    
    global nDataPoints;
    nDataPoints = nDataPoints+batchsize;
end

