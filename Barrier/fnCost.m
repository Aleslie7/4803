function  [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x, u, k,R,dt, Q, p_target)
u = u - 1.2263;
l0 = 0.5 * u' * R *u + 0.5 * (x - p_target)' * Q *(x - p_target);
l_x = Q * (x - p_target);
l_xx = Q;
l_u = R * u;
l_uu = R;
l_ux = zeros(1,13);

end