dx = 0.1;
l = 3;
m1 = 1.5;
m2 = 1;
e = 0.4*dx*dx;
T = 4;
dt = 0.002;
N = 2*l/dx;

syms x y;
W12(x,y) = (x^2+y^2)/2;
W21(x,y) = -W12(x,y);
W11(x,y) = (x^2+y^2)/2;
W22(x,y) = W11(x,y);
W.W11 = W11; W.W22 = W22; W.W21 = W21; W.W12 = W12;

X = -l+dx/2:dx:l-dx/2;
X = meshgrid(X);
Y = X.';

r1_0 = zeros(N,N);
r1_0(-1<X & X<0 & -0.5<Y & Y<0.5) = m1;
r2_0 = zeros(N,N);
r2_0(0<X & X<1 & -0.5<Y & Y<0.5) = m2;

r1 = r1_0;
r2 = r2_0;

[r1,r2] = two2d ( r1,r2,l,W ,dt, T, e,'v');
