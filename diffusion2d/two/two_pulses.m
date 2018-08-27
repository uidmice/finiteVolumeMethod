dx = 0.1;
l = 3;
m1 = 1;
m2 = 5;
% e = 0;
e = 0.4*dx*dx;
nt = 2000;
dt = 0.002;
N = 2*l/dx;

syms x y H;
W12(x,y) = sqrt(x^2+y^2);
W21(x,y) = -W12(x,y);
W11(x,y) = (x^2+y^2)/2;
W22(x,y) = W11(x,y);
H(x) = 0;
V = zeros(N,N);



X = -l+dx/2:dx:l-dx/2;
X = meshgrid(X);
Y = X.';

r1_0 = zeros(N,N);
r1_0(-1<X & X<0 & -0.5<Y & Y<0.5) = m1;
r2_0 = zeros(N,N);
r2_0(0<X & X<1 & -0.5<Y & Y<0.5) = m2;

r1 = r1_0;
r2 = r2_0;
r = zeros(N,N,2);
for i = 1:nt
    [r1,r2] = two2d (l, r1,r2, W11,W12,W22,W21, H, V, dt, dt, e);
    if any(r1<0 | r2<0)
        return
    end
    s = sprintf('attr_attr%03d.mat',i);
    r(:,:,1) = r1;
    r(:,:,2) = r2;
    save(s,'r');
end