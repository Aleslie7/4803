global Horizon;
global time;
global p_target;
global dt;


%% Plot Convergence Information.
   figure(1);
   subplot(3,2,1)
   hold on
   plot(time(1:end-1),x_traj(1,1:end-1),'linewidth',4);  
   %plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
   title('X Velocity','fontsize',20); 
    xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   
    subplot(3,2,2);hold on;
   plot(time(1:end-1),x_traj(2,1:end-1),'linewidth',4); 
   %plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
   title('Y Velocity','fontsize',20);
    xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
    
    subplot(3,2,3);hold on
   plot(time(1:end-1),x_traj(3,1:end-1),'linewidth',4); 
   %plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)
   title('Omega','fontsize',20)
     xlabel('Time in sec','fontsize',20)
   hold off;
   grid;


%    subplot(3,2,4);hold on
%    plot(x_traj(1,:),x_traj(2,:),'linewidth',2); 
%    title('Vehicle Trajectory','fontsize',20)
%      xlabel('X Coordinate','fontsize',20)
%      ylabel('Y Coordinate')
%    hold off;
%    grid;
% Cost(3,13:100) = 0;
   
    subplot(3,2,5);hold on
   plot(Cost(3,:),'linewidth',2); 
    xlabel('Iterations','fontsize',20)
   title('Cost','fontsize',20);
   %save('DDP_Data');
   time2 = 1:1:149;
    if (length(u_k) < 100)
        time2 = 1:3:147
    end
   subplot(3,2,6);hold on
   plot(time2, u_k,'linewidth',2); 
    xlabel('Input','fontsize',20)
   title('Controller','fontsize',20);





