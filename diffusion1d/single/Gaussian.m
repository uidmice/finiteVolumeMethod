dx = 0.1;
l = 20;
x = -l+dx/2:dx:l-dx/2;
r0 = zeros(1,length(x));
r0(x<1/2 & x>-1/2)=exp(-1./(1-x(x<1/2 & x>-1/2).^2))/2;
T = 200;
dt = 0.001;
syms W1(y) W2(y) W3(y)
W1(y) = exp(-y^2/2); %||W1||=2.5
W2(y) = exp(-y^2/2)/2;
W3(y) = exp(-y^2/2)*5;
ita1 = 1;
ita2 = 2;
[r1, E1] = single1d (r0,l, W1, dt, T, ita1);
[r2, E2] = single1d (r0,l, W1, dt, T, ita2);

         