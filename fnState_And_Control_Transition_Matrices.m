function [A,B] = fnState_And_Control_Transition_Matrices(x,u,du,dt, dynamics)
x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);

u1 = u(1,1);
u2 = u(2,1);

A = zeros(3,3);
B = zeros(3,2);
A(1,3) = -u1.*sin(x3);
A(2,3) = u1.*cos(x3);
B(1,1) = cos(x3);
B(2,1) = sin(x3);
end