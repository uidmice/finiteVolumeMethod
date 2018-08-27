dx = 0.1;
l = 3;
m1 = 0.6;
m2 = 0.1;
e = 0.4*dx^2;
nt = 200;
dt = 0.002;
N = 2*l/dx;

syms x y H;
W12(x,y) = sqrt(x^2+y^2);
W21(x,y) = W12(x,y);
W11(x,y) = (x^2+y^2)/2;
W22(x,y) = W11(x,y);
H(x) = 0;
V = zeros(N,N);



X = -l+dx/2:dx:l-dx/2;
X = meshgrid(X);
Y = X.';

dis = sqrt(X.^2+Y.^2);
r1_0 = zeros(N,N);
r1_0(dis<2) = m1/(pi*4);
r2_0 = zeros(N,N);
r2_0(dis<0.5) = m2/(pi*0.5*0.5);

r1 = r1_0;
r2 = r2_0;
r = zeros(N,N,2);
subplot(1,2,1)
surf(X,Y,r1_0);
subplot(1,2,2)
surf(X,Y,r2_0);
drawnow
for i = 1:nt
    [r1,r2] = two2d (l, r1,r2, W11,W12,W22,W21, H, V, dt, dt, e);
    s = sprintf('attr_attr%03d.mat',i);
    r(:,:,1) = r1;
    r(:,:,2) = r2;
    save(s,'r');
    subplot(1,2,1)
    surf(X,Y,r1);
    subplot(1,2,2)
    surf(X,Y,r2);
    title(sprintf('t=%d', i));
    drawnow
end