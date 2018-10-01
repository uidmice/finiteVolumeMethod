function [dr1,dr2, E1,E2] = drho(r1,r2,K11,K12,K21,K22,V,dx,  ita)

r = [r1;r2];
ks1 = -(dx*conv(r1,K11,'same')-dx*conv(r2,K12,'same')+ V);
ks2 = -(dx*conv(r2,K22,'same')-dx*conv(r1,K21,'same')+ V);

u = -(circshift(ks1,-1)-ks1)/dx;
v = -(circshift(ks2,-1)-ks2)/dx;

theta = 2;
rhox = px(r,dx,theta);
rhoE = r+dx*rhox/2;
rhoW = r-dx*rhox/2;

Fx = max(u,0).*rhoE(1,:)+min(u,0).*circshift(rhoW(1,:),-1);
Gx = max(v,0).*rhoE(2,:)+min(v,0).*circshift(rhoW(2,:),-1);

dr1 = -(Fx-circshift(Fx,1))/dx + ita * lap(r1.*r1, dx);
dr2 = -(Gx-circshift(Gx,1))/dx + ita * lap(r2.*r2, dx);
E1 = sum(ks1.*r1)*dx/2;
I1 = sum(r1.*r1)*ita*dx;
E1 = E1 + I1;

E2 = sum(ks2.*r2)*dx/2;
I2 = sum(r2.*r2)*ita*dx;
E2 = E2 + I2;
end


function rhox = px(rho, dx, theta)
rho_p = circshift(rho,-1);
rho_m = circshift(rho,1);
r1 = theta*(rho-rho_m);
r2 = (rho_p-rho_m)/2;
r3 = theta*(rho_p-rho);
rhox = arrayfun(@minmod,r1,r2,r3)/dx;
end

function r = lap(x, dx)
r = (circshift(x,1)-2*x+circshift(x,-1))/dx^2;
end

