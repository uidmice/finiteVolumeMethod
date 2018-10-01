%Define 
l = 4;                                  
dx = 0.1;
dt = 0.01;
T = 2;
X = meshgrid(-l+dx/2:dx:l-dx/2);
Y = X.';
N = length(X);

rho0 = zeros(N,N);
rho0(X<3 & X>-3 & Y<3 & Y>-3) = 1/4;
syms W(x,y) H(x)
W(x,y) = -exp(- x^2 - y^2)/pi;
H(x) = x^3/30;

rhoR = single2d (rho0,l, W, dt, T, 'v', H);