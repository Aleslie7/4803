function [dynamics] = fnDynamics()
syms x [12 1]
syms u [4 1]

% static variables
m = 0.5; % kg
Ixx = 0.0032; % kgm^2
Iyy = 0.0032; % kgm^2
Izz = 0.0055; % kgm^2
kt = 0.01691; % 1/m
l = 0.17; % m

% doing this here for easy debugging
R = [cos(x8)*cos(x9) cos(x8)*cos(x9) -sin(x8);
    sin(x8)*sin(x7)*cos(x9)-cos(x7)*sin(x9) sin(x8)*sin(x7)*cos(x9)+cos(x7)*cos(x9) sin(x7)*cos(x8);
    sin(x7)*sin(x9)+cos(x7)*sin(x8)*cos(x9) sin(x8)*sin(x9)*cos(x7)-sin(x7)*sin(x9) cos(x8)*cos(x7)];
n = [0; 0; (u1 + u2 + u3 + u4)/m];
temp = R*n;

A = [0 0 0 1 0 0 0 0 0 0 0 0; % x_dot
    0 0 0 0 1 0 0 0 0 0 0 0; % y_dot
    0 0 0 0 0 1 0 0 0 0 0 0; % z_dot
    0 0 0 0 0 0 0 0 0 0 0 0; % x_ddot
    0 0 0 0 0 0 0 0 0 0 0 0; % y_ddot
    0 0 0 0 0 0 0 0 0 0 0 0; % z_ddot
    0 0 0 0 0 0 0 0 0 1 tan(x8)*sin(x7) tan(x8)*cos(x7); % phi_dot
    0 0 0 0 0 0 0 0 0 0 cos(x7) -sin(x7); % theta_dot
    0 0 0 0 0 0 0 0 0 0 sin(x7)/cos(x8) cos(x7)/cos(x8); % psi_dot
    0 0 0 0 0 0 0 0 0 0 0 -(Izz - Iyy)*x11/Ixx; % p_dot
    0 0 0 0 0 0 0 0 0 0 0 (Izz - Ixx)*x10/Iyy; % q_dot
    0 0 0 0 0 0 0 0 0 0 0 0]; % r_dot

B = [0 0 0 0; % x_dot
    0 0 0 0; % y_dot
    0 0 0 0; % z_dot
    0 0 0 0; % x_ddot these arent actually 0 but inside of temp
    0 0 0 0; % y_ddot
    0 0 0 0; % z_ddot
    0 0 0 0; % phi_dot
    0 0 0 0; % theta_dot
    0 0 0 0; % psi_dot
    0.7071*l/Ixx -0.7071*l/Ixx 0.7071*l/Ixx -0.7071*l/Ixx; % p_dot
    -0.7071*l/Iyy -0.7071*l/Iyy 0.7071*l/Iyy 0.7071*l/Iyy; % q_dot
    kt/Izz -kt/Izz -kt/Izz kt/Izz]; % r_dot

% adding rotation matrix stuff in
Fa = A * x(1:12) + [0; 0; 0; 0; 0; -9.81; 0; 0; 0; 0; 0; 0];
G = B * u + [0; 0; 0; temp(1); temp(2); temp(3); 0; 0; 0; 0; 0; 0];
F = Fa + G;


% creating symbolic stuff for ease of interaction later
dynamics = struct();
dynamics.F = matlabFunction(F);
dynamics.Fx = matlabFunction(jacobian(F, x));
dynamics.Fu = matlabFunction(jacobian(F, u));
dynamics.Fb =  matlabFunction(G);
end