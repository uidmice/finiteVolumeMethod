function [rhop, E] = drho(rho,K,V,dx,  ita)

ks = 2*ita*rho-(dx*conv(rho,K,'same')+ V);

u = -(circshift(ks,-1)-ks)/dx;

theta = 2;
rhox = px(rho,dx,theta);
rhoE = rho+dx*rhox/2;
rhoW = rho-dx*rhox/2;

Fx = max(u,0).*rhoE+min(u,0).*circshift(rhoW,-1);

rhop = -(Fx-circshift(Fx,1))/dx ;
E = sum(ks.*rho)*dx/2;
end


function rhox = px(rho, dx, theta)
rho_p = circshift(rho,-1);
rho_m = circshift(rho,1);
r1 = theta*(rho-rho_m);
r2 = (rho_p-rho_m)/2;
r3 = theta*(rho_p-rho);
rhox = arrayfun(@minmod,r1,r2,r3)/dx;
end

