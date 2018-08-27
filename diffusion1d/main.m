
m1 = 0.6;
m2 = 0.1;
e = 1.7;
N = 80;
buff = 10;


[x, rho0, eta0, W11,W12,W21,W22] = initialize("AA",3,m1,m2,e,1.5,buff,N);
xprime = x(1:end-1)+(x(2)-x(1))/2;