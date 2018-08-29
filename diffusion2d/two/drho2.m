function [dr] = drho2(r, K11,K12,K21,K22, dH, V, dx, e)

r1 = r(:,:,1);
r2 = r(:,:,2);
pHij = double(dH(r1+r2));

ks1 = dx*dx*(conv2(r1,K11,'same')+conv2(r2,K12,'same'))+ pHij+V+e*(r1+r2);
ks2 = dx*dx*(conv2(r2,K22,'same')+conv2(r1,K21,'same'))+ pHij+V+e*(r1+r2);

ks = cat(3,ks1,ks2);
u = -(circshift(ks,-1)-ks)/dx;
v = -(circshift(ks,-1,2)-ks)/dx;

theta = 2;
rhox = px(r,dx,theta);
rhoy = py(r,dx,theta);
rhoE = r + dx*rhox/2;
rhoW = r - dx*rhox/2;
rhoN = r + dx*rhoy/2;
rhoS = r - dx*rhoy/2;

Fx1 = max(u(:,:,1),0).*rhoE(:,:,1)+min(u(:,:,1),0).*circshift(rhoW(:,:,1),-1);
Fx2 = max(u(:,:,2),0).*rhoE(:,:,2)+min(u(:,:,2),0).*circshift(rhoW(:,:,2),-1);
Fy1 = max(v(:,:,1),0).*rhoN(:,:,1)+min(v(:,:,1),0).*circshift(rhoS(:,:,1),-1,2);
Fy2 = max(v(:,:,2),0).*rhoN(:,:,2)+min(v(:,:,2),0).*circshift(rhoS(:,:,2),-1,2);

Fx = cat(3,Fx1,Fx2);
Fy = cat(3,Fy1,Fy2);

dr = -(Fx-circshift(Fx,1)+Fy-circshift(Fy,1,2))/dx;

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

