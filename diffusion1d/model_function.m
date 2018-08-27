
function [rhoR, etaR] = model_function (x, rho0, eta0, W11,W12,W21,W22, e)
N = length(x)-1;
dx = (x(N+1)-x(1))/N;
[U,V] = getDensityVelocity(rho0,eta0,dx,W11,W12,W21,W22,e);

persistent i r
if isempty(i)
    i = 1;
    r = ones(1,500);
end


r(i) = dx/(max(max(U),max(V))*4);
if i>=500
    i = 1;
else 
    i = i + 1;
end
dt = min(r);

[F,G] = getFlux (U,V,rho0,eta0);
[rhoR, etaR] = step (rho0, eta0, F, G, dt, dx);

end




function [U,V] = getDensityVelocity (rho,eta,dx,W11,W12,W21,W22, e)

%{
input:
    rho,eta:        densities   1*N
    dx:             step size
    W:              interaction potential functions
    e:              coefficient of the cross-diffusivity    
output:
    U,V:            density velocity 1*(N+1)
%}
    N= length(rho);
    t1 = zeros(1,N+1);
    t2 = zeros(1,N+1);
    
%     ind = linspace(0,1-N,N);
    
    
%     for j = 1:N
%         t1(j)=(dot(rho,W11(ind*dx))+dot(eta,W12(ind*dx)))*dx+e*(rho(j)+eta(j));
%         t2(j)=(dot(eta,W22(ind*dx))+dot(rho,W21(ind*dx)))*dx+e*(rho(j)+eta(j));
%         ind = ind + 1;
%     end
    
    mat = zeros(N,N);
    for i = 1:N
        mat(:,i) = -(1-i:N-i)';
    end
    
    mat = mat*dx;
    t1(1:N) = dx*(rho*double(W11(mat))+eta*double(W12(mat)))+e*(rho+eta);
    t2(1:N) = dx*(eta*double(W22(mat))+rho*double(W21(mat)))+e*(rho+eta);
    t1(end) = t1(1);
    t2(end) = t2(1);
    
    U = zeros(1,N+1);
    V = zeros(1,N+1);
    
    U(2:end) = (t1(1:N)-t1(2:N+1))/dx;
    V(2:end) = (t2(1:N)-t2(2:N+1))/dx;
    U(1)=U(end);
    V(1)=V(end);
    
end


function [F,G] = getFlux (U,V,rho,eta)
%{
input:
    rho,eta:        densities           1*N
    U,V:            density velocity    1*(N+1)
output:
    F,G:    flux .   1*(N+1)
%}

F = max(0,U)*diag([0,rho])+min(0,U)*diag([rho,0]);
G = max(0,V)*diag([0,eta])+min(0,V)*diag([eta,0]);

end


function [rhoN, etaN] = step (rho, eta, F, G, dt, dx)
%{
input:
    rho,eta:    densities       1*N
    F,G:        flux .          1*(N+1)
    dt:         step size
    dx:         space size
output:
    rhoN,etaN:  densities at next time step     1*N
%}

rhoN = rho - (F(2:end)-F(1:end-1))*dt/dx;
etaN = eta - (G(2:end)-G(1:end-1))*dt/dx;

end





