
function [rhoR] = single2d (varargin)
%SINGLE2D: numeric simulation of the dynamics of one single species in 2D
%   r = single2d(r0,l, W, dt, T) 
%   Input:
%       r0:         initial density, NxN matrix
%       l:          domain[-l,l]x[-l,l]
%       W:          W(x) interacting potentials
%       dt:         time step
%       T:          simulation time. Total #iterations = T/dt
%   Output: 
%       r:          density at t = T, NxN matrix
%
%   Optional parm:
%   r = single2d(.. ,H): H is a symbolic function for internal energy as a
%   function of the density. Default H(r) = 0
%
%   r = single2d(.. ,V) optionally sets the environmental confinement
%   potential V, which is a NxN matrix. Default: 0
%
%   r = single2d(.. ,e) sets the diffusion coefficient for some e > 0.
%   Default e = 0
%
%   r = single2d(.. ,'v') or r = single2d(.. ,'V') enables visual display
%   during the simulation. Default disabled.
%   
%   r = single2d(.. ,'solver') where 'solver' sets the numeric method used
%   for ODE. Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.
%   

narginchk(5,10);


[rho0,l, W, dt, T, e, V, H ,v, s]  = arguments(varargin{:});

N = length(rho0);
dx = 2*l/N;
K = getKernal(W, dx, N);
dH = diff(H);
L = @(X) drho(X,K,dH,V,dx, e);
nt = ceil(T/dt);
rhoR = rho0;
X = -l+dx/2:dx:l-dx/2;
X = meshgrid(X);
Y = X.';
for i = 1:nt
    if s=='s'
        rhoR = SSPRK3(rhoR,L,dt);
    else
        rhoR = Euler(rhoR,L,dt);
    end 
    
    if v
        draw(X,Y,rhoR);
        title(sprintf('t = %.3f', i*dt));
        drawnow;
    end
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

function [rho0,l, W, dt, T, e, V, H ,v, s]   = arguments(varargin)

rho0 = varargin{1};
l = varargin{2};
W = varargin{3};
dt = varargin{4};
T = varargin{5};

N = length(rho0);
syms x H(x)

% Optional input defaults
e = 0;                      % No diffusion term
v = false;                  % No visible update during simulation
V = zeros(N, N);            % No confinement potential
H(x) = 0;                   % No internal energy
s = 's';

% Loop over optional arguments
for i = 6:nargin
    a = varargin{i};
    if ischar(a) && isscalar(a) && lower(a) == 'v'
        v = true;
    elseif isreal(a) && isscalar(a) && isfinite(a)
        e = a;
    elseif isa(a,'symfun')
        H = a;
    elseif isstring(a) && ~strcmp(a,'Euler')
        s = 'e';
    elseif ismatrix(a) && length(a)==N
        V = a;
    elseif ismatrix(a) && length(a)~=N
        error("The shape of the input confinement energy matrix is incompatible");
    else
        error('arguments:nonsense','Failed to interpret argument #%d.',i);
    end
                
end

end






