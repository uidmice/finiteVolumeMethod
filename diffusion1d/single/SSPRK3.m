function [u, E] = SSPRK3(u0,L,dt)
[u1, E] = L(u0);
u1 = u0+dt*u1;
u2 = 3/4*u0+u1/4+dt*L(u1)/4;
u = u0/3+2/3*u2+2/3*dt*L(u2);
end

