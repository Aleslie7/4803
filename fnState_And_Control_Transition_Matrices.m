function [A,B] = fnState_And_Control_Transition_Matrices(x,u,du,dt, dynamics)
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

A = dynamics.Fx(u1, u2, u3, u4, x7, x8, x10, x11, x12);
B = dynamics.Fu(x7, x8);
end