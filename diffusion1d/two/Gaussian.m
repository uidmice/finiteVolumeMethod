dx = 0.1;
l = 10;
x = -l+dx/2:dx:l-dx/2;
r1 = zeros(1,length(x));
r1(x<1/2 & x>-1/2)=exp(-1./(1-x(x<1/2 & x>-1/2).^2))/2;
r2 = r1.*2;
T = 200;
dt = 0.001;
syms W11(y) W12(y) W22(y) W21(y)
W11(y) =  exp(-y^2/2); %||W1||=2.5
W22(y) = W11(y);
W12(y) = exp(-y^2/2)/2;
W21(y) = W12(y);
W.W11 = W11; W.W22 = W22; W.W21 = W21; W.W12 = W12;

ita = 4;
[R1, R2, E1,E2] = two1d (r1,r2,l,W,dt,T,ita);
