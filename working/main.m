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
Q(1,1)=1; %X
Q(2,2)=1; %Y
Q(3,3)=120; %Z
Q(4,4)=.3; %Xdot
Q(5,5)=.3; %Ydot
Q(6,6)=.3; %Zdot
Q(7,7)=1; %Phi
Q(8,8)=1; %Theta
Q(9,9)=1; %Psy
Q(10,10)=50; %P
Q(11,11)=50; %q
Q(12,12)=50; %r

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

% Learning Rate
gamma = 0.3;
dynamics = fnDynamics();

for k = 1:num_iter
%------------------------------------------------> Linearization of the dynamics
%------------------------------------------------> Quadratic Approximations of the cost function 
    for  j = 1:(Horizon-1)    
        [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j), j,R,dt, Q, p_target);
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
    Cost(:,k) = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R, Q);
    
    if (k ~= 1)
        residuals(:, k) = abs(Cost(:, k - 1) - Cost(:, k));
    end
%     x1(k,:) = x_traj(1,:);

    if mod(k, 10) == 0
        fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
    end

end

residuals(:,1) = residuals(:,2);
x_traj(:,end) - p_target

time = zeros(1,Horizon);
time(1)=0;
for i= 2:Horizon
    time(i) =time(i-1) + dt;  
end
