
function [r1,r2] = two2d (l, r1_0,r2_0, W11,W12,W22,W21, H, V, dt, T, e)
%input:
%   l:          domain[-l,l]x[-l,l]
%   r1_0,r2_0:  initial density
%   W's:        W(x) interacting potentials
%   H:          H(rho) density of internal energy
%   V:          V(x) confinement potential
%   dt:         step size
%   T:          simulation time
%   e:          diffusion coef

N = length(r1_0);
dx = 2*l/N;
K11 = getKernal(W11, dx, N);
K12 = getKernal(W12, dx, N);
K21 = getKernal(W21, dx, N);
K22 = getKernal(W22, dx, N);
dH = diff(H);
L = @(X) drho2(X, K11,K12,K21,K22, dH, V, dx, e);
nt = ceil(T/dt);
R = zeros([N,N,2]);
R(:,:,1) = r1_0;
R(:,:,2) = r2_0;

for i = 1:nt
    R = SSPRK3(R,L,dt);
end

r1 = R(:,:,1);
r2 = R(:,:,2);

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







