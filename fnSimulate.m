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
    u1 = u_new(1,k);
    u2 = u_new(2,k);
    u3 = u_new(3,k);
    u4 = u_new(4,k);

    x(:,k+1) = x(:,k) + (dynamics.F(u1, u2, u3, u4, x4, x5, x6, x7, x8, x10, x11, x12) + [0; 0; 0; 0; 0; -9.81; 0; 0; 0; 0; 0; 0])* dt; % adds gravity after as static value
end
