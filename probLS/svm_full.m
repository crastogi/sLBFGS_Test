function [f, df, vf2, vdf2, var_f, var_df, grad_matrix] = svm_full(x, sampleidx, varargin)
    global data classes batchsize nDataPoints;
    lambda = .001;
    if ~exist('sampleidx','var')
        sampleidx = randsample(1:(size(data, 1)), batchsize);
    end
    n = length(sampleidx);
   
    % x is col vector
    % -- function value ---------------------------------------------------
    % vector-compute alpha value
    alpha = sum(data(sampleidx,:).*x', 2);
    % Vector-compute temp value
    temp = max(0, 1-classes(sampleidx).*alpha);
    % Vector-compute function value
    f = sum(temp.^2)/(2*n) + lambda*(x'*x)/2;
    %lf = var(temp.^2/(2));
    
    % Mass-compute gradient
%    classes(sampleidx).*alpha < 1 -> then compute gradient
    compute_gradient = (classes(sampleidx).*alpha<1)*1; % gives logical array
    % y*(1-y)
    grad_scalar = classes(sampleidx).*(1-classes(sampleidx).*alpha);
    % Set to zero the values that we do not want to update
    grad_scalar = grad_scalar.*compute_gradient;
    grad_matrix = data(sampleidx,:);
    % Row-wise multiply by scalar values (-1*y*(1-y))
    for i = 1:n
        grad_matrix(i,:) = -1*grad_matrix(i,:)*grad_scalar(i) + lambda*x';
    end
    % Compute overall gradient, normalize by batch size and include the
    % regularizer; Note that output needs to be a column vector
    df = sum(grad_matrix,1)/n;
    df = df';
    
    % Build vf and vdf
    var_f = 1/(n-1)*(sum(temp.^4/4)/n - (sum(temp.^2/2)/n)^2);
    var_df = 1/(n-1)*(sum(grad_matrix.^2,1)/n - (sum(grad_matrix,1)/n).^2);
    % output needs to be a column vector
    var_df = var_df';
    
    vf2 = var_f;
    vdf2 = var_df;

    nDataPoints = nDataPoints + n;
    
    return;
end