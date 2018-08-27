
function [rhoR] = single2d (l, rho0, W, H, V, dt, T, e)
%input:
%   l:          domain[-l,l]x[-l,l]
%   rho0:       initial density
%   W:          W(x) interacting potentials
%   H:          H(rho) density of internal energy
%   V:          V(x) confinement potential
%   dt:         step size
%   T:          simulation time
%   e:          diffusion coef
N = length(rho0);
dx = 2*l/N;
K = getKernal(W, dx, N);
dH = diff(H);
L = @(X) drho(X,K,dH,V,dx, e);
nt = ceil(T/dt);
rhoR = rho0;
for i = 1:nt
    rhoR = SSPRK3(rhoR,L,dt);
end

end

function [K] = getKernal (W, dx, N)

%{
input:
    W:              interaction potential function
    dx:             space size
    N:              number of points 
output:
    Wij:            help matrix for calculation W(i, j)
%}

   xi = 0:N-1;
   xi = repmat(xi, 1,N)*dx;
   yi = 0:N-1;
   yi = repelem(yi, N)*dx;
   
   K = double(W(xi, yi));
   
   t = logical(isnan(K)+isinf(K));
   if any(t)
       in = find(double(t));
       a = 0.5773503*dx/2;
       if length(in)==1
           K(in) =(W(xi(in)+a,yi(in)+a)+W(xi(in)-a,yi(in)+a)+W(xi(in)+a,yi(in)-a)+W(xi(in)-a,yi(in)-a))/4;
       else
           for i = range (length(in))
               ind = in(i);
               K(ind) =(W(xi(ind)+a,yi(ind)+a)+W(xi(ind)-a,yi(ind)+a)+W(xi(ind)+a,yi(ind)-a)+W(xi(ind)-a,yi(ind)-a))/4;
           end
       end
   end
   K = reshape(K, N, N);
   
   hf = flip(K(:,2:end),2);
   vf = flip(K(2:end,:));
   f = flip(vf(:,2:end),2);
   K = [f, vf;hf,K];
end







