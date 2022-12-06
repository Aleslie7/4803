function [dynamics] = fnDynamics()
x = sym('x', [1 3]);
u = sym('u',[1 2]);

Fa = [0;0;0];

G = [u(1).*cos(x(3));
    u(1).*sin(x(3));
    u(2)];

F = Fa + G;

dynamics = struct();
dynamics.F = matlabFunction(F);
dynamics.Fx = matlabFunction(jacobian(F, x));
dynamics.Fu = matlabFunction(jacobian(F, u));
dynamics.Fb = matlabFunction(G);
end
