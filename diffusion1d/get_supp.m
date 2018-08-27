function [supp_r,supp_e] = get_supp(m1,m2,e)
%Input:
%       m1, m2:         mass (for rho and eta respectively)
%       e:              coefficient of the cross-diffusivity                
%Output:
%       supp_r, supp_eta:   support range (c, and b)

syms b c;
eqn1 = 3*e*(m1+m2*b)==3*e*(m2+m1*b)+(c-b)^2*(3*m2+2*c*m1+m1*b);
eqn2 = m1*(c^2-b^2)+e*(m1+m2)+2*m2*(c-b)==2*sqrt(e)*(m2+m1*b)*cot(b/sqrt(e));
range = [0 10 ; 0 10];

supp_r = 0;
supp_e = 1;
while supp_r<supp_e
    [sol_b, sol_c] = vpasolve([eqn1, eqn2], [b, c], range, 'random',true);
    supp_r = double(sol_c);
    supp_e = double(sol_b);
end
    
end

