function [u, cost] = fnDDP(x,num_iter, dt, Q_f, R, p_target, gamma,...
    ~, x_dim, u_dim, u_init, dynamics,Q)
global Horizon;

xo = x;
x_traj = zeros(x_dim,Horizon); 
if nargin < 11
 u_k = zeros(u_dim,Horizon-1);
else
 u_k = u_init;
end

du_k = zeros(u_dim,Horizon-1);
    
for k = 1:num_iter
    q0 = zeros(Horizon-1);
    q_k = zeros(x_dim, Horizon-1);
    Q_k = zeros(x_dim, x_dim, Horizon-1);
    r_k = zeros(u_dim, Horizon-1);
    R_k = zeros(u_dim, u_dim, Horizon-1);
    P_k = zeros(u_dim, x_dim, Horizon-1);
    A = zeros(x_dim, x_dim, Horizon-1);
    B = zeros(x_dim, u_dim, Horizon-1);
    for  j = 1:(Horizon-1)    
        [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j),...
            j,R,dt,Q,p_target);
        q0(j) = dt * l0;        % L.
        q_k(:,j) = dt * l_x;    % Lx.
        Q_k(:,:,j) = dt * l_xx; % Lxx.
        r_k(:,j) = dt * l_u;    % Lu.
        R_k(:,:,j) = dt * l_uu; % Luu.
        P_k(:,:,j) = dt * l_ux; % Lux.
        [fx,fu] = fnState_And_Control_Transition_Matrices(x_traj(:,j),...
            u_k(:,j),du_k(:,j),dt, dynamics);
        A(:,:,j) = eye(x_dim,x_dim) + fx * dt;
        B(:,:,j) = fu * dt;  
    end

    Vxx = zeros(x_dim,x_dim,Horizon);
    Vx = zeros(x_dim, Horizon);
    V = zeros(1, Horizon);
    Vxx(:,:,Horizon)= Q_f;
    Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target); 
    V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f *...
        (x_traj(:,Horizon) - p_target);

    for j = (Horizon-1):-1:1
        H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j);
        G = P_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j);
        g_ = r_k(:,j) +  B(:,:,j)' * Vx(:,j+1);

        inv_H = inv(H); 
        L_k(:,:,j)= - inv_H * G; 
        l_k (:,j) = - inv_H *g_;  
        Vxx(:,:,j) = Q_k(:,:,j)+ A(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) +...
            L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j);
        Vx(:,j)= q_k(:,j) +  A(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g_...
            + G' * l_k(:,j) + L_k(:,:,j)'*H * l_k(:,j);
        V(:,j) = q0(j) + V(j+1) + 0.5 *  l_k (:,j)' * H * l_k (:,j) +...
            l_k (:,j)' * g_;
    end
    u_new = zeros(u_dim, Horizon-1);
    dx = zeros(x_dim,1);

    for i=1:(Horizon-1)
        du = l_k(:,i) + L_k(:,:,i) * dx;
        dx = A(:,:,i) * dx + B(:,:,i) * du;  
        u_new(:,i) = u_k(:,i) + gamma * du;
    end

    u_k = u_new;
    [x_traj] = fnSimulate(xo,u_new,Horizon,dt,0, dynamics);
end
cost = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R,Q);
u = u_k;
end