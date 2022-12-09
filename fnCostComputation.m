function [Cost] =  fnCostComputation(x_traj,u_new,p_target,dt,Q_f,R, Q)
    [numOfStates,Horizon] = size(x_traj);
    Cost = 0;

    for j =1:(Horizon-1)

        %Cost = Cost + (0.5 * x_traj(:,j)' * Q * x_traj(:,j) + 0.5 * u_new(:,j)' * R * u_new(:,j)) * dt;
        Cost = Cost + (0.5 * (x_traj(:,j) - p_target)' * Q * (x_traj(:,j) - p_target) + 0.5 * u_new(:,j)' * R * u_new(:,j)) * dt;

    end

    TerminalCost= (x_traj(:,Horizon) - p_target)'*Q_f * (x_traj(:,Horizon) - p_target); 
    Cost = Cost + TerminalCost;
   