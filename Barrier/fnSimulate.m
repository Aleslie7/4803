function [x] = fnSimulate(xo,u_new,Horizon,dt,sigma, dynamics, num_timesteps)

% checks if normal DDP or MPC step
steps = Horizon - 1;
if exist('num_timesteps') == 1
    steps = num_timesteps;
end

x = xo;

for k = 1:steps
    x1 = x(1,k);
    x2 = x(2,k);
    x3 = x(3,k);
    x4 = x(4,k);
    x5 = x(5,k);
    x6 = x(6,k);
    x7 = x(7,k);
    x8 = x(8,k);
    x9 = x(9,k);
    x10 = x(10,k);
    x11 = x(11,k);
    x12 = x(12,k);
    wk = x(13,k);
    u1 = u_new(1,k);
    u2 = u_new(2,k);
    u3 = u_new(3,k);
    u4 = u_new(4,k);

        % barrier states for beta_k
    h1 = (x1 - 2.2)^2 + (x2 - 2.2)^2 + (x3 - 1)^2 - 1;
    h2 = (x1)^2 + (x2 - 0.2)^2 + (x3)^2 - 1;
    h3 = (x1 - 3)^2 + (x2)^2 + (x3 - 0.5)^2 - 1;
    bk = 1/h1 + 1/h2 + 1/h3;
    
    % using desired states for beta_d
    h1 = (5 - 2.2)^2 + (3 - 2.2)^2 + (2 - 1)^2 - 1;
    h2 = (5)^2 + (3 - 0.2)^2 + (2)^2 - 1;
    h3 = (5 - 3)^2 + (3)^2 + (2 - 0.5)^2 - 1;
    bd = 1/h1 + 1/h2 + 1/h3;
    
    wk = bk - bd; % w_k

    temp = dynamics.F(u1, u2, u3, u4, x1, x2, x3, x4, x5, x6, x7, x8, x10, x11, x12);
    x(1:12,k+1) = x(1:12,k) + (temp(1:12))* dt; % 1:12 to exclude stepping of w_k

    % updating for beta_k+1
    x1 = x(1,k+1);
    x2 = x(2,k+1);
    x3 = x(3,k+1);
    h1 = (x1 - 2.2)^2 + (x2 - 2.2)^2 + (x3 - 1)^2 - 1;
    h2 = (x1)^2 + (x2 - 0.2)^2 + (x3)^2 - 1;
    h3 = (x1 - 3)^2 + (x2)^2 + (x3 - 0.5)^2 - 1;
    bk1 = 1/h1 + 1/h2 + 1/h3;

    gamma = 0.5; % usually 0.3-0.7
    wk1 = bk1 - bd - gamma * (wk + bd - bk); % w_k+1, gamma term will 0 out since we have no noise
    x(13,k+1) = wk1;


end