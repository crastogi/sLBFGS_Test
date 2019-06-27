% Follows algorithm 1
% M: max epochs
% N: number of data points
% d: number of dimensions
% k: `k' in k-SVRG
% eta: step size
% Phi: tracks selected indices
% aBar_m: 'overall' gradient
k = 5;
bs = 25;
M = 1000;
eta = .25;
epsilon = 1E-5;
global nDataPoints
nDataPoints = 0;

N   = size(data,1);
d   = size(data,2);
l   = ceil(N/(k*bs));
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
for m = 0:(M-1)                             % Line 4
    Phi = [];                               % Line 5
	x_t_pos = zeros(l,d);
    for t = 0:(l-1)                         % Line 6
        x_t_pos(t+1,:) = x_t;           
        i_t = randsample(1:N, bs);          % Line 7
        % Line 8
        sub_bindings = theta_m_binding(i_t);
        a_it = 0;
        for curr_theta = unique(sub_bindings)
            [~, sg] = func(theta_m(curr_theta,:)', i_t(sub_bindings==curr_theta));
            a_it = a_it + sum(sub_bindings==curr_theta)*sg;
        end
        a_it = a_it/bs;
        [~, f_it] = func(x_t', i_t);
        g_it = f_it' - a_it' + aBar_m';
        x_t = x_t - eta*g_it;               % Line 9
        Phi = [Phi, i_t];          % Line 10
    end
    Phi = unique(Phi);
    xBar = sum(x_t_pos, 1)/S_l;             % Line 12
    % New theta_i: Line 17; start by setting theta_m_binding
    theta_m1_binding = theta_m_binding;
    theta_m1_binding(Phi) = max(theta_m1_binding)+1;
    min_Tidx = min(theta_m1_binding);
    max_Tidx = max(theta_m1_binding);
    % Reallocate the theta matrix
    theta_m1 = zeros((max_Tidx-min_Tidx+1), d);
    for curr_Tidx = min_Tidx:(max_Tidx-1)
        theta_m1(curr_Tidx-min_Tidx+1,:) = theta_m(curr_Tidx,:);
    end
    theta_m1(max_Tidx-min_Tidx+1,:) = xBar;
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
    theta_m_binding = theta_m1_binding-min(theta_m1_binding)+1;
    theta_m = theta_m1;
    disp([num2str(log10(norm(aBar_m))) ' ' num2str(log10(norm(xBar'-x_min{1}))) ' ' num2str(sum(theta_m_binding==1)*(sum(theta_m(1,:))==0))]);
    if (norm(aBar_m)/max([1 norm(xBar)]) < epsilon)
        disp("Convergence!");
        return;
    end
end