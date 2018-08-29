l = 2;
dx = 0.1;
dt = 0.01;
T = 2;
X = meshgrid(-l+dx/2:dx:l-dx/2);
Y = X.';
N = length(X);

rho0 = zeros(N,N);
rho0(X<1 & X>-1 & Y<1 & Y>-1) = 1;
syms W(x,y)
W(x,y) = x^2/2 - log(x^2 + y^2)/2 + y^2/2;

rhoR = single2d (rho0,l, W, dt, T, 'v', 0); %set diffusion coeff to 0

rhoRp = single2d (rho0,l, W, dt, T, 'v', 0.4*dx^2); %set diffusion coeff to 0.4*dx^2