function  [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x, u, k,R,dt)

l0 = u' *R *u;
l_x =zeros(3,1);
l_xx = zeros(3,3);
l_u = R * u;
l_uu = R;
l_ux = zeros(2,3);

end