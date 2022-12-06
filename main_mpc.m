clear;
close all;

global Horizon;
global time;
global p_target
global dt;

% Environment parameters.


% Obtain expressions for F, Fx, Fu.
dynamics_nominal = fnDynamics();

% Environment parameters with measurement error.
dynamics_sigma = 0.2;


%dynamics_actual = fnDynamics(mc_noisy, mp_noisy, l_noisy, g_noisy);
dynamics_actual = fnDynamics();

% Solver parameters.
Horizon = 150;  % Time Horizon.
num_iter = 100; % Number of Iterations
dt = 0.02;     % Discretization.

% Costs.
Q_f = zeros(3,3); % State cost. 4x4 since state is 4-dimensional.
Q_f(1,1) = 10;     % X position cost.
Q_f(2,2) = 50;   % X velocity cost.
Q_f(3,3) = 200;  % Pole angle cost.


if ~(all(eig(Q_f) >= 0))
    error('Cost matrix Q_f not positive semi-definite.')
end

R = 1* eye(2,2); % Control cost. 1x1 since control is 1-dimensional.

% Initialize solution.
% State represented as [x, x_dot, theta, theta_dot].
xo = [-1.5;2;0];
x_dim = length(xo);
u_dim = 2;

% Goal state:
p_target = zeros(x_dim, 1);
p_target(1,1) = 2.5; % Target x
p_target(2,1) = -1.5; % Target x_dot
p_target(3,1) = pi./2; % Target theta


% Add noise.
sigma_nominal = 0.0;
sigma_real = 0.0;

% Learning Rate
gamma = 0.5;

% Run the MPC.
x = xo;
residuals = []; % Residual history.
Cost = [];      % Cost history.
x_traj = [];    % Initial trajectory
u_init = zeros(u_dim, Horizon-1);

i = 1;
max_num_iters = 150;
while 1
    [u_k, cost] = fnDDP(x,num_iter, dt, Q_f, R, p_target, gamma,...
      ~, x_dim, u_dim, u_init, dynamics_nominal);
    u_init = u_k;
    Cost = [Cost cost];
    [x_new] = fnSimulate(x,u_k,Horizon,dt,sigma_real, dynamics_actual, 2);
    x = x_new(:,2);
    x_traj = [x_traj x];
    if i ~= 1
        residuals = [residuals abs(Cost(:, i - 1) - Cost(:, i))];
    end
    if mod(i, 100) == 0
        fprintf('MPC Iteration %i,  Current Cost = %e \n',i,...
            norm(x - p_target));
    end
    
    if (norm(x - p_target) < 1e-3) || (i > max_num_iters)
        break;
    end
    i = i + 1;
end

Horizon = length(x_traj);
time = zeros(1, Horizon);
for i= 2:(Horizon - 1)
    time(i) =time(i-1) + dt;
end

visualize