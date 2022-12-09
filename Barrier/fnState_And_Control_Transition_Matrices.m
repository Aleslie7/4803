function [A,B] = fnState_And_Control_Transition_Matrices(x,u_new,du,dt, dynamics)
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
wk = x(13,1);
u1 = u_new(1,1);
u2 = u_new(2,1);
u3 = u_new(3,1);
u4 = u_new(4,1);

A = dynamics.Fx(u1, u2, u3, u4, x1, x2, x3, x7, x8, x10, x11, x12);
B = dynamics.Fu(x7, x8);
end