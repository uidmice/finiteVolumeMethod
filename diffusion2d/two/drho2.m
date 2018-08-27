function [rhop] = drho2(rho, K11,K12,K21,K22, dH, V, dx, e)

r1 = rho(:,:,1);
r2 = rho(:,:,2);
pHij = double(dH(r1+r2));

ks1 = dx*dx*(conv2(r1,K11,'same')+conv2(r2,K12,'same'))+ pHij+V+e*(r1+r2);
ks2 = dx*dx*(conv2(r2,K22,'same')+conv2(r1,K21,'same'))+ pHij+V+e*(r1+r2);

ks = cat(3,ks1,ks2);
u = -(circshift(ks,-1)-ks)/dx;
v = -(circshift(ks,-1,2)-ks)/dx;
% u1 = -(circshift(ks1,-1)-ks1)/dx;
% u2 = -(circshift(ks2,-1)-ks2)/dx;
% v1 = -(circshift(ks1,-1,2)-ks1)/dx;
% v2 = -(circshift(ks2,-1,2)-ks2)/dx;

theta = 2;
rhox = px(rho,dx,theta);
rhoy = py(rho,dx,theta);
rhoE = rho + dx*rhox/2;
rhoW = rho - dx*rhox/2;
rhoN = rho + dx*rhoy/2;
rhoS = rho - dx*rhoy/2;
% rhox1 = px(r1,dx,theta);
% rhox2 = px(r2,dx,theta);
% rhoy1 = py(r1,dx,theta);
% rhoy2 = py(r2,dx,theta);
% rhoE1 = r1+dx*rhox1/2;
% rhoE2 = r2+dx*rhox2/2;
% rhoW1 = r1-dx*rhox1/2;
% rhoW2 = r2-dx*rhox2/2;
% rhoN1 = r1+dx*rhoy1/2;
% rhoN2 = r2+dx*rhoy2/2;
% rhoS1 = r1-dx*rhoy1/2;
% rhoS2 = r2-dx*rhoy2/2;

Fx1 = max(u(:,:,1),0).*rhoE(:,:,1)+min(u(:,:,1),0).*circshift(rhoW(:,:,1),-1);
Fx2 = max(u(:,:,2),0).*rhoE(:,:,2)+min(u(:,:,2),0).*circshift(rhoW(:,:,2),-1);
Fy1 = max(v(:,:,1),0).*rhoN(:,:,1)+min(v(:,:,1),0).*circshift(rhoS(:,:,1),-1,2);
Fy2 = max(v(:,:,2),0).*rhoN(:,:,2)+min(v(:,:,2),0).*circshift(rhoS(:,:,2),-1,2);

Fx = cat(3,Fx1,Fx2);
Fy = cat(3,Fy1,Fy2);

rhop = -(Fx-circshift(Fx,1)+Fy-circshift(Fy,1,2))/dx;
% rhop(:,:,1) = -(Fx1-circshift(Fx1,1)+Fy1-circshift(Fy1,1,2))/dx;
% rhop(:,:,2) = -(Fx2-circshift(Fx2,1)+Fy2-circshift(Fy2,1,2))/dx;

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

