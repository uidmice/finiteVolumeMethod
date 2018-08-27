function [u] = Euler(u0,L,dt)
u = u0+dt*L(u0);
end

