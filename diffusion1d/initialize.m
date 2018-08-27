function [x, rho0, eta0, W11,W12,W21,W22,rho_a, eta_a] = initialize(type,cond,m1,m2,e, l,buff, N)
%Input:
%       type:   model           'AA/AR'("attractive-attractive/attractive-repulsive")
%       cond:   conditions      
%                       1: Batman profiles (symmetric, stationary)
%       m1, m2: mass
%       e:      cross-diffusion coeffi
%       l:      x domain [-l,l]
%       buff:   number of empty buffer points each side around [-l,l]
%       N:      number of space divisions on [-l,l]
%Output: 
%       x:              support points                                      1*(N+1+2*buff)
%       rho0,eta0:      initial densities (based on 'type' and 'cond') .    1*(N+2*buff)
%       W:              interaction potentials
%       rho_a, eta_a:   analysic solutions (optional)
syms y
W11(y) = y^2/2;
W22(y) = y^2/2;
W12(y) = abs(y);
if type=="AA"
    W21(y) = abs(y);
else
    W21(y) = -abs(y);
end

dx = 2 * l / N;
L = l + buff*dx;
x = linspace(-L,L,N+1+2*buff);
xm = x(1:end-1)+dx/2;

rho0 = zeros(1,N+2*buff);
eta0 = zeros(1,N+2*buff);
switch cond
    case 1
        rho0(xm>=-l & xm<=l) = m1/(2*l);
        eta0(xm>=-l*0.8 & xm<=l*0.8) = eta0 + m2/(2*l*0.8);
        [c, b] = get_supp(m1,m2,e);
        u = (m2+m1*b)/(sqrt(e)*sin(b/sqrt(e)));
        
        f1(y) = -m1*(y^2-c^2)/(2*e)+m2*(y+c)/e;
        f2(y) = u*cos(y/sqrt(e))/2-m2/2;
        f3(y) = -m1*(y^2-c^2)/(2*e)-m2*(y-c)/e;
        rho_a(y) = piecewise(-c<=y<-b,f1(y), -b<=y<=b, f2(y), b<y & y<=c, f3(y),0 );
        eta_a(y) = piecewise(-b<=y<=b,u*cos(y/sqrt(e))/2-m1/2, 0);        
    case 3
        rho0(abs(xm)<=l/3) = m1/(2/3*l);
        eta0((abs(xm)>=l/3 & abs(xm)<=l))= m2/(4/3*l);
end





end

