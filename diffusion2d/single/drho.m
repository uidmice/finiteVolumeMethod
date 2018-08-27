function [rhop] = drho(rho, K, dH, V, dx, e)

pHij = double(dH(rho));
ks = dx*dx*conv2(rho,K,'same')+ pHij+V+e*rho;

u = -(circshift(ks,-1)-ks)/dx;
v = -(circshift(ks,-1,2)-ks)/dx;

theta = 2;
rhox = px(rho,dx,theta);
rhoy = py(rho,dx,theta);
rhoE = rho+dx*rhox/2;
rhoW = rho-dx*rhox/2;
rhoN = rho+dx*rhoy/2;
rhoS = rho-dx*rhoy/2;

Fx = max(u,0).*rhoE+min(u,0).*circshift(rhoW,-1);
Fy = max(v,0).*rhoN+min(v,0).*circshift(rhoS,-1,2);

rhop = -(Fx-circshift(Fx,1)+Fy-circshift(Fy,1,2))/dx;

end


function rhox = px(rho, dx, theta)
rho_up = circshift(rho,-1);
rho_down = circshift(rho,1);
r1 = theta*(rho-rho_down);
r2 = (rho_up-rho_down)/2;
r3 = theta*(rho_up-rho);
rhox = arrayfun(@minmod,r1,r2,r3)/dx;
end

function rhoy = py(rho, dx, theta)
rho_right = circshift(rho,1,2);
rho_left = circshift(rho,-1,2);
r1 = theta*(rho-rho_right);
r2 = (rho_left-rho_right)/2;
r3 = theta*(rho_left-rho);
rhoy = arrayfun(@minmod,r1,r2,r3)/dx;
end

