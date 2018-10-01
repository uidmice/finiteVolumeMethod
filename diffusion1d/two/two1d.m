
function [R1, R2, E1,E2] = two1d (varargin)
%TWO2D: numeric simulation of the dynamics of two species in 1D
%   [R1, R2, E1,E2] = two2d(r1, r2, l, W, dt, T) 
%   Input:
%       r1, r2:     initial density, 1xN arrays
%       l:          domain[-l,l]
%       W:          W11(x), W12(x), W21(x), W22(x) interacting potentials
%       dt:         time step
%       T:          simulation time. Total #iterations = T/dt
%   Output: 
%       R1,R2:      densities from t =0 to T, (T/dt) x N matrix
%       E1,E2:      energies from t = 0 to T, 1 x (T/dt) arrays
%
%   Optional parm:
%
%   [R1, R2, E1,E2] = single1d(.. ,V) optionally sets the environmental confinement
%   potential e * V(x), which is a 1xN array. Default: 0
%
%   [R1, R2, E1,E2] = single1d(.. ,ita) sets the diffusion coeffi. Default ita = 1
%
%   [R1, R2, E1,E2] = single1d(.. ,'v') or r = single2d(.. ,'V') enables visual display
%   during the simulation. Default disabled.
%   
%   [R1, R2, E1,E2] = single1d(.. ,'solver') where 'solver' sets the numeric method used
%   for ODE. Possible options: 'Euler', 'SSPRK3'. Default 'SSPRK3'.
%   

narginchk(6,10);


[r1,r2, l, W, dt, T, ita, V,v, s]  = arguments(varargin{:});

N = length(r1);
dx = 2*l/N;
K11 = getKernal(W.W11, dx, N);
K12 = getKernal(W.W12, dx, N);
K21 = getKernal(W.W21, dx, N);
K22 = getKernal(W.W22, dx, N);

L = @(X1,X2) drho(X1,X2,K11,K12,K21,K22,V,dx,  ita);
nt = ceil(T/dt);
R1 = zeros(nt+1,N);
R2 = zeros(nt+1,N);

E1 = zeros(1,nt+1);
E2 = zeros(1,nt+1);
[R1(1,:),R2(1,:), E1(1),E2(1)] = L(r1,r2);
R1(1,:) = r1;
R2(1,:) = r2;
for i = 1:nt
    if s=='s'
        [R1(i+1,:),R2(i+1,:), E1(i+1),E2(i+1)]= SSPRK3(R1(i,:),R2(i,:),L,dt);
    else
        [R1(i+1,:), E1(i+1)] = Euler(R1(i,:),L,dt);
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

function [r1,r2, l, W, dt, T, ita, V,v, s] = arguments(varargin)

r1 = varargin{1};
r2 = varargin{2};
l = varargin{3};
W = varargin{4};
dt = varargin{5};
T = varargin{6};

N = length(r1);
syms x 

% Optional input defaults
ita = 1;                      % No diffusion term
v = false;                  % No visible update during simulation
V = zeros(1, N);            % No confinement potential                  % No internal energy
s = 's';

% Loop over optional arguments
for i = 7:nargin
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






