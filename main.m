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
Q=eye(12,12);
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

% Weight in Final State:
Q_f = 10*eye(12,12);
Q_f(1,1)=100;
Q_f(2,2)=100;
Q_f(3,3)=150;

% Weight in the Control:
R = diag([0.001 0.001 0.001 0.001]);

% Initial configuration:
xo = [-3; -3; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
x_dim = length(xo);

% Initial Trajectory:
x_traj = zeros(x_dim,Horizon);

% Initial Control:
%u_k = zeros(2,Horizon-1);
u_k = zeros(4,Horizon-1);
u_dim = size(u_k, 1);
du_k = zeros(4,Horizon-1);
%du_k = zeros(u_dim,Horizon-1);

%Cost = zeros(3, num_iter); % Cost history.
%residuals = zeros(3, num_iter); % Residual history.

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

% Learning Rate
gamma = 0.5;
dynamics = fnDynamics();

for k = 1:num_iter
    % Preallocate cost memory.
    q0 = zeros(Horizon-1); % (1,)
    q_k = zeros(x_dim, Horizon-1);
    Q_k = zeros(x_dim, x_dim, Horizon-1);
    r_k = zeros(u_dim, Horizon-1);
    R_k = zeros(u_dim, u_dim, Horizon-1);
    P_k = zeros(u_dim, x_dim, Horizon-1);
    A = zeros(x_dim, x_dim, Horizon-1);
    B = zeros(x_dim, u_dim, Horizon-1);
    for  j = 1:(Horizon-1)    
        [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j), j,R,dt);
        % Compute loss function gradients for the current timestep.
        % Quadratic Approximations of the cost function.
        q0(j) = dt * l0; % L.
        q_k(:,j) = dt * l_x; % Lx.
        Q_k(:,:,j) = dt * l_xx; % Lxx.
        r_k(:,j) = dt * l_u; % Lu.
        R_k(:,:,j) = dt * l_uu; % Luu.
        P_k(:,:,j) = dt * l_ux; % Lux.

        % Linearize the dynamics using first order taylor expansion.
        [fx,fu] = fnState_And_Control_Transition_Matrices(x_traj(:,j),u_k(:,j),du_k(:,j),dt, dynamics);
        A(:,:,j) = eye(x_dim,x_dim) + fx * dt;
        B(:,:,j) = fu * dt;  
    end

    % Preallocate value function memory.
    Vxx = zeros(x_dim,x_dim,Horizon);
    Vx = zeros(x_dim, Horizon);
    V = zeros(1, Horizon);
 
    % Compute value function at final timestep, its gradient, and its jacobian.
    Vxx(:,:,Horizon)= Q_f;
    Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target); 
    V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target);
    
    % Backpropagation of the Value Function.
    for j = (Horizon-1):-1:1
        % Quu = Luu + B^T * Vxx * B
        H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j);
        % Qux = Lux + B^T * Vxx * A
        G = P_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j);
        % Qu = Lu + B^T * Vx
        g_ = r_k(:,j) +  B(:,:,j)' * Vx(:,j+1);

        thing = eye(length(H));
        inv_H = H \ thing; % Quu^-1
        L_k(:,:,j)= - inv_H * G; % Feedback term = -Quu^-1 * Qux.
        l_k (:,j) = - inv_H *g_; % Feedforward term = -Quu^-1 * Qu.

        % Vxx = (Lxx + A^T * Vxx * A) + (-Qxu * Quu^-1) * Quu * (-Quu^-1 * Qux)
        % + (-Qxu * Quu^-1 * Qux) + (Qxu * -Quu^-1 * Qux) 
        Vxx(:,:,j) = Q_k(:,:,j)+ A(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) + L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j);
        % Vx = (Lx + A^T * Vx') + (-Qxu * Quu^-1 * Qu) + (Qxu * -Quu^-1 * Qu) +
        % (-Qxu * Quu^-1 * Qu * Quu * -Quu^-1 * Qu)
        Vx(:,j)= q_k(:,j) +  A(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g_ + G' * l_k(:,j) + L_k(:,:,j)'*H * l_k(:,j);
        % V = L + V' + (0.5 * Qu^T * Quu^-1 * Quu * Quu^-1 * Qu) + (-Qu^T * Quu^-1 * Qu)
        V(:,j) = q0(j) + V(j+1) + 0.5 *  l_k (:,j)' * H * l_k (:,j) + l_k (:,j)' * g_;
    end

    % Preallocate control variation memory.
    u_new = zeros(u_dim, Horizon-1);
    dx = zeros(x_dim,1);
%     for i=1:(Horizon)
for i=1:(Horizon-1)
        % Find the controls.
        du = l_k(:,i) + L_k(:,:,i) * dx;
        dx = A(:,:,i) * dx + B(:,:,i) * du;  
        u_new(:,i) = u_k(:,i) + gamma * du;
        
        % prevent control from going negative
        u_new(:,i) = max(u_new(:,i), [0;0;0;0]);
    end

    u_k = u_new;




    [x_traj] = fnSimulate(xo,u_new,Horizon,dt,0, dynamics); %change above variables for new dynamics to test controller
    Cost(:,k) = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
    %     x_trag 3x150
    %     u_k 2x149
    %     p_target 3x1
    %     dt 1x1
    %     Q_f 3x3
    
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

% if ~visualizing_bundles
%     visualize;
% %     close(fh);
% end