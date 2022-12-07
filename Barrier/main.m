visualizing_bundles_exists = exist('visualizing_bundles', 'var');
visualizing_bundles = visualizing_bundles_exists && visualizing_bundles;

if ~visualizing_bundles
    close all;
    visualizing_bundles=false;
end

global Horizon;
global time;
global p_target
global Q_f;
global dt;

% Horizon 
Horizon = 800; % 1.5sec
% Number of Iterations
num_iter = 100;

% Discretization
dt = 0.01;

% Weight running:
Q=eye(13,13);
Q(1,1)=.1; %X
Q(2,2)=.1; %Y
Q(3,3)=120; %Z
Q(4,4)=.1; %Xdot
Q(5,5)=.1; %Ydot
Q(6,6)=.1; %Zdot
Q(7,7)=1; %Phi
Q(8,8)=1; %Theta
Q(9,9)=1; %Psy
Q(10,10)=70; %P
Q(11,11)=70; %q
Q(12,12)=70; %r
Q(13,13)=10; %barrier

% Weight in Final State:
Q_f = 10*eye(13,13);
Q_f(1,1)=100;
Q_f(2,2)=100;
Q_f(3,3)=150;

% Weight in the Control:
R = diag([0.001 0.001 0.001 0.001]);

% Initial barrier state finding:
% barrier states for beta_k
h1 = (-3 - 2.2)^2 + (-3 - 2.2)^2 + (-1 - 1)^2 - 1;
h2 = (-3)^2 + (-3 - 0.2)^2 + (-1)^2 - 1;
h3 = (-3 - 3)^2 + (-3)^2 + (-1 - 0.5)^2 - 1;
bk = 1/h1 + 1/h2 + 1/h3;

% using desired states for beta_d
h1 = (5 - 2.2)^2 + (3 - 2.2)^2 + (2 - 1)^2 - 1;
h2 = (5)^2 + (3 - 0.2)^2 + (2)^2 - 1;
h3 = (5 - 3)^2 + (3)^2 + (2 - 0.5)^2 - 1;
bd = 1/h1 + 1/h2 + 1/h3;

wk = bk - bd; % w_k


% Initial configuration:
xo = [-3; -3; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; wk];
x_dim = length(xo);

% Initial Trajectory:
x_traj = zeros(x_dim,Horizon);

% Initial Control:
%u_k = zeros(2,Horizon-1);
u_k = zeros(1,Horizon-1);
u_dim = size(u_k, 1);
du_k = zeros(1,Horizon-1);
%du_k = zeros(u_dim,Horizon-1);

%Cost = zeros(3, num_iter); % Cost history.
%residuals = zeros(3, num_iter); % Residual history.

% Final barrier state finding:
% Initial barrier state finding:
% barrier states for beta_k
h1 = (5 - 2.2)^2 + (3 - 2.2)^2 + (2 - 1)^2 - 1;
h2 = (5)^2 + (3 - 0.2)^2 + (2)^2 - 1;
h3 = (5 - 3)^2 + (3)^2 + (2 - 0.5)^2 - 1;
bk = 1/h1 + 1/h2 + 1/h3;

wk = bk - bd; % w_k


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
p_target(12,1) = wk;

% Learning Rate
gamma = 0.5;
dynamics = fnDynamics();

for k = 1:num_iter
%------------------------------------------------> Linearization of the dynamics
%------------------------------------------------> Quadratic Approximations of the cost function 
    for  j = 1:(Horizon-1)    
        [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j), j,R,dt);
        q0(j) = dt * l0; % L.
        q_k(:,j) = dt * l_x; % Lx.
        Q_k(:,:,j) = dt * l_xx; % Lxx.
        r_k(:,j) = dt * l_u; % Lu.
        R_k(:,:,j) = dt * l_uu; % Luu.
        P_k(:,:,j) = dt * l_ux; % Lux.

        [dfx,dfu] = fnState_And_Control_Transition_Matrices(x_traj(:,j),u_k(:,j),du_k(:,j),dt, dynamics);

        A(:,:,j) = eye(x_dim,x_dim) + dfx * dt;
        B(:,:,j) = dfu * dt;  
    end
 
%------------------------------------------------> Find the controls
    Vxx(:,:,Horizon)= Q_f;
    Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target); 
    V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target);
    
%------------------------------------------------> Backpropagation of the Value Function
    for j = (Horizon-1):-1:1
        H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j);
        G = P_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j);
        g_d = r_k(:,j) +  B(:,:,j)' * Vx(:,j+1);

        thing = eye(length(H));
        inv_H = H \ thing; % Quu^-1
        L_k(:,:,j)= - inv_H * G; % Feedback
        l_k (:,j) = - inv_H *g_d; % Feedforward

        Vxx(:,:,j) = Q_k(:,:,j)+ A(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) + L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j);
        Vx(:,j)= q_k(:,j) +  A(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g_d + G' * l_k(:,j) + L_k(:,:,j)'*H * l_k(:,j);
        V(:,j) = q0(j) + V(j+1) + 0.5 *  l_k (:,j)' * H * l_k (:,j) + l_k (:,j)' * g_d;
    end

%----------------------------------------------------> Forward Propagaion: 

%--------------------------------> Find the controls/ forward 
    dx = zeros(x_dim,1);
for i=1:(Horizon-1)
        % Find the controls.
        du = l_k(:,i) + L_k(:,:,i) * dx;
        dx = A(:,:,i) * dx + B(:,:,i) * du;  
        u_new(:,i) = u_k(:,i) + gamma * du;

        % prevent control from going negative
        u_new(:,i) = max(u_new(:,i), [0;0;0;0]);
    end

    u_k = u_new;



%-------------------------------> Simulation of the Nonlinear System
    [x_traj] = fnSimulate(xo,u_new,Horizon,dt,0, dynamics); %change above variables for new dynamics to test controller
    Cost(:,k) = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
    
    if (k ~= 1)
        residuals(:, k) = abs(Cost(:, k - 1) - Cost(:, k));
    end
%     x1(k,:) = x_traj(1,:);

    if mod(k, 10) == 0
        fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
    end

end

residuals(:,1) = residuals(:,2);

time = zeros(1,Horizon);
time(1)=0;
for i= 2:Horizon
    time(i) =time(i-1) + dt;  
end

if ~visualizing_bundles
    visualize;
%     close(fh);
end
