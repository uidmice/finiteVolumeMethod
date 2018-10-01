
function [rhoR, E] = single1d (varargin)
%SINGLE2D: numeric simulation of the dynamics of one single species in 2D
%   r = single2d(r0,l, W, dt, T) 
%   Input:
%       r0:         initial density, 1xN array
%       l:          domain[-l,l]
%       W:          W(x) interacting potentials
%       dt:         time step
%       T:          simulation time. Total #iterations = T/dt
%   Output: 
%       r:          density at t = T, 1xN array
%
%   Optional parm:
%
%   r = single1d(.. ,V) optionally sets the environmental confinement
%   potential e * V(x), which is a 1xN array. Default: 0
%
%   r = single1d(.. ,ita) sets the diffusion coeffi. Default ita = 1
%
%   r = single1d(.. ,'v') or r = single2d(.. ,'V') enables visual display
%   during the simulation. Default disabled.
%   
%   r = single1d(.. ,'solver') where 'solver' sets the numeric method used
%   for ODE. Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.
%   

narginchk(5,9);


[r0,l, W, dt, T, ita, V,v, s]  = arguments(varargin{:});

N = length(r0);
dx = 2*l/N;
K = getKernal(W, dx, N);
L = @(X) drho(X,K,V,dx,  ita);
nt = ceil(T/dt);
rhoR = zeros(nt+1,N);
E = zeros(1,nt+1);
[rhoR(1,:), E(1)] = drho(r0,K,V,dx,  ita);
rhoR(1,:) = r0;
for i = 1:nt
    if s=='s'
        [rhoR(i+1,:), E(i+1)]= SSPRK3(rhoR(i,:),L,dt);
    else
        [rhoR(i+1,:), E(i+1)] = Euler(rhoR(i,:),L,dt);
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
    K:              kernal
%}

   xi = 0:N-1;
   K = double(W(xi));
   
   t = logical(isnan(K)+isinf(K));
   if any(t)
       in = find(double(t));
       a = 0.1*dx/2;
       if length(in)==1
           K(in) =(W(xi(in)+a)+W(xi(in)-a))/2;
       else
           for i = range (length(in))
               ind = in(i);
               K(ind) =(W(xi(ind)+a)+W(xi(ind)-a))/2;
           end
       end
   end
   
   f = flip(K(2:end));
   K = [f,K];
end

function [rho0,l, W, dt, T, ita, V ,v, s]   = arguments(varargin)

rho0 = varargin{1};
l = varargin{2};
W = varargin{3};
dt = varargin{4};
T = varargin{5};

N = length(rho0);
syms x 

% Optional input defaults
ita = 1;                      % No diffusion term
v = false;                  % No visible update during simulation
V = zeros(1, N);            % No confinement potential                  % No internal energy
s = 's';

% Loop over optional arguments
for i = 6:nargin
    a = varargin{i};
    if ischar(a) && isscalar(a) && lower(a) == 'v'
        v = true;
    elseif isreal(a) && isscalar(a) && isfinite(a)
        ita = a;
    elseif isstring(a) && ~strcmp(a,'Euler')
        s = 'e';
    elseif isvector(a) && length(a)==N
        V = a;
    elseif isvector(a) && length(a)~=N
        error("The size of the input confinement energy matrix is incompatible");
    else
        error('arguments:nonsense','Failed to interpret argument #%d.',i);
    end
                
end

end






