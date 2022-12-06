function [A,B] = fnState_And_Control_Transition_Matrices(x,u,du,dt, dynamics)
x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);
x5 = x(5,1);
x6 = x(6,1);
x7 = x(7,1);
x8 = x(8,1);
x9 = x(9,1);
x10 = x(10,1);
x11 = x(11,1);
x12 = x(12,1);

u1 = u(1,1);
u2 = u(2,1);
u3 = u(3,1);
u4 = u(4,1);

A = zeros(12,12);
B = zeros(12,4);

% Needs to be updated
A(1,3) = -u1.*sin(x3);
A(2,3) = u1.*cos(x3);
B(1,1) = cos(x3);
B(2,1) = sin(x3);
end
