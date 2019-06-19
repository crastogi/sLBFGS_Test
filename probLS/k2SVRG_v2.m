% Follows algorithm 2
% M: max epochs
% N: number of data points
% d: number of dimensions
% k: `k' in k-SVRG
% eta: step size
% Phi: tracks selected indices
% aBar_m: 'overall' gradient
% -------NOTE: Works pretty well??
k = 5;
M = 1000;
eta = .1;
epsilon = 1E-5;
global nDataPoints;
nDataPoints = 0;

N   = size(data,1);
d   = size(data,2);
l   = ceil(N/k);
x0  = flip(x_min{1}')*5; %zeros(1,d);
S_l = l;                                    % Line 3 
% Need an size N array to store idx of locations. Initialize to 1
theta_m_binding = ones(1, N);
% Initialize theta matrix to 1 position set to x0
theta_m = x0;
% Set aBar_m
[~, aBar_m] = func(x0', 1:N);
x_t = x0;

% Outer loop
for m = 0:(M-1)                             % Line 5
    ind = randperm(N);                      % Line 6
    theta_m1_binding = theta_m_binding;
    theta_m1 = theta_m;
    for j = 0:(k-1)                         % Line 7
%        Phi = [];                           % Line 8
        x_t_pos = zeros(l,d);
        for t = 0:(l-1)                     % Line 9
            x_t_pos(t+1,:) = x_t;
            i_t = randi(N);                 % Line 10
            % Line 11
            [~, a_it] = func(theta_m(theta_m_binding(i_t),:)', i_t);
            [~, f_it] = func(x_t', i_t);
            g_it = f_it' - a_it' + aBar_m';
            x_t = x_t - eta*g_it;           % Line 12
        end
        Phi = ind((j*l+1):(j*l+l));         % Line 13
        xBar = sum(x_t_pos, 1)/S_l;         % Line 15
        % New theta_i: Line 17; start by setting theta_m_binding
        theta_m1_binding(Phi) = max(theta_m1_binding)+1;
        % Reallocate the theta matrix
        theta_m1(max(theta_m1_binding),:) = xBar;
    end
    % Re-number the theta_m1 matrix
    min_Tidx = min(theta_m1_binding);
    max_Tidx = max(theta_m1_binding);
    theta_temp = zeros((max_Tidx-min_Tidx+1), d);
    for curr_Tidx = min_Tidx:max_Tidx
        theta_temp(curr_Tidx-min_Tidx+1,:) = theta_m1(curr_Tidx,:);
    end
    % Need to update full gradient estimate: Line 18
    [~, new_phi_grad] = func(xBar', Phi);
    % Compute update components from previous theta_m's
    old_bindings = theta_m_binding(Phi);
    old_grad = 0;
    for curr_theta = unique(old_bindings)
        [~, og] = func(theta_m(curr_theta,:)', Phi(old_bindings==curr_theta));
        old_grad = old_grad + sum(old_bindings==curr_theta)*og;
    end
    aBar_m = aBar_m + (new_phi_grad*length(Phi) - old_grad)/N;
    % Now adjust theta_m_binding indices
    theta_m_binding = theta_m1_binding-min_Tidx+1;
    theta_m = theta_temp;
    if (norm(aBar_m)/max([1 norm(xBar)]) < epsilon)
        disp("Convergence!");
        return;
    end
    disp([num2str(m) ' ' num2str(log10(norm(aBar_m))) ' ' num2str(log10(norm(xBar'-x_min{1}))) ' ' num2str(sum(theta_m_binding==1)*(sum(theta_m(1,:))==0))]);
end