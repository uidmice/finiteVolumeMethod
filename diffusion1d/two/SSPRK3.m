function [U1,U2, E1,E2] = SSPRK3(u1,u2,L,dt)
[U1,U2, E1,E2] = L(u1,u2);
U1 = u1+dt*U1;
U2 = u2+dt*U2;
[K1,K2] = L(U1,U2);
K1 = 3/4*u1+U1/4+dt*K1/4;
K2 = 3/4*u2+U2/4+dt*K2/4;

[P1,P2] = L(K1,K2);
U1 = u1/3+2/3*K1+2/3*dt*P1;
U2 = u2/3+2/3*K2+2/3*dt*P2;
end

