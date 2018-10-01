function [u, E] = Euler(u0,L,dt)
[u, E] = L(u0);
u = u0+dt*u;
end

