clear;
close all;

global Horizon;
global time;
global p_target
global dt;

% Environment parameters.
Q=eye(12,12);
Q(1,1)=1; %X
Q(2,2)=1; %Y
Q(3,3)=120; %Z
Q(4,4)=.3; %Xdot
Q(5,5)=.3; %Ydot
Q(6,6)=20; %Zdot
Q(7,7)=1; %Phi
Q(8,8)=1; %Theta
Q(9,9)=1; %Psy
Q(10,10)=50; %P
Q(11,11)=50; %q
Q(12,12)=50; %r

% Obtain expressions for F, Fx, Fu.
dynamics_nominal = fnDynamics();

% Environment parameters with measurement error.
dynamics_sigma = 0;


%dynamics_actual = fnDynamics(mc_noisy, mp_noisy, l_noisy, g_noisy);
dynamics_actual = fnDynamics();

% Solver parameters.
Horizon = 75;  % Time Horizon.
num_iter = 10; % Number of Iterations
dt = 0.01;     % Discretization.

% Costs.
Q_f = 10*eye(12,12);
Q_f(1,1)=100;
Q_f(2,2)=100;
Q_f(3,3)=150;


if ~(all(eig(Q_f) >= 0))
    error('Cost matrix Q_f not positive semi-definite.')
end

R = diag([0.001 0.001 0.001 0.001]);

% Initialize solution.

% Initial configuration:
xo = [-3; -3; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
x_dim = length(xo);

% Initial Trajectory:
x_traj = zeros(x_dim,Horizon);

% Initial Control:
u_k = zeros(4,Horizon-1) + 1.2263;
u_dim = size(u_k, 1);
du_k = zeros(4,Horizon-1);

% Target:
p_target(1,1) = 5;
p_target(2,1) = 3;
p_target(3,1) = 2;
p_target(4,1) = 0;
p_target(5,1) = 0;
p_target(6,1) = 0;
p_target(7,1) = 0;
p_target(8,1) = 0;
p_target(9,1) = 0;
p_target(10,1) = 0;
p_target(11,1) = 0;
p_target(12,1) = 0;
% Add noise.
sigma_nominal = 0.0;
sigma_real = 0.0;

% Learning Rate
gamma = 0.3;

% Run the MPC.
x = xo;
residuals = []; % Residual history.
Cost = [];      % Cost history.
x_traj = [];    % Initial trajectory
u_init = zeros(u_dim, Horizon-1);

i = 1;
max_num_iters = 30;
while 1
    [u_k, cost] = fnDDP(x,num_iter, dt, Q_f, R, p_target, gamma,...
      0, x_dim, u_dim, u_init, dynamics_nominal,Q);
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